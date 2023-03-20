library(R2jags)
library(mvtnorm) # for multivariate normal distribution
library(ggplot2)
library(patchwork)

## -------------------------------------------------------------------------- ##
## ------------------ Multivariate model (env) --------------------- ##
## -------------------------------------------------------------------------- ##

rm(list = ls())

genData <- function(N, G, I, NE, muT, sigma, sg, se, st){

  g <- rnorm(I, 0, sg)
  e <- rnorm(G, 0, se)
  gen <- rep(1:I, each = N/I)

  # cluster membership
  env <- rep(NA, N)
  t <- matrix(NA, ncol = NE, nrow = N)
  for (i in 1:N) {
    env[i] <- sample(1:G, size = 1, replace = TRUE)
    t[i, ] <- rmvnorm(1, muT[env[i], ], st)
  }
  # env <- sample(1:G, size = N, replace = TRUE)
  # t <- rnorm(N, muT[env], st)

  # generate mean of the model
  mu_y <- g[gen] + e[env]
  y <- rnorm(N, mu_y, sigma)

  df <- data.frame(y = y, env = env, gen = gen)

  return(list(df = df, g = g, e = e, t = t))

}


G <- 3
N <- 100
NE <- 2
muT <- matrix(c(-5, 0, 2, -5, 0, 2), ncol = 2, nrow = G)
sigma <- 1
sg <- 10
se <- 10
st <- diag(2)
I <- 5

set.seed(02)
dat <- genData(N = N, G = G, NE = NE, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st)
plot(dat$t)


# Jags code to fit the model to the simulated data
model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] #env here is a parameter but gen is not

    # # Clustering model
    t[i, 1:NE] ~ dmnorm(muT[env[i], 1:NE], st_inv)
    env[i] ~ dcat(pi[i, 1:G])

    for (g in 1:G) {
      exp_theta[i, g] <- exp(theta[i, g])
      pi[i, g] <- exp(theta[i, g]) / sum(exp_theta[i, 1:G])
      theta[i, g] ~ dnorm(0, 6^-2)
    }

   }

  # # Prior on cluster membership
  # for (g in 1:G) {
  #     pi[1:NE, g] ~ ddirch(alpha)
  #   }


  # Priors on means
  for (g in 1:G) {
    for (m in 1:NE) {
      muT_raw[g, m] ~ dnorm(0, 100^-2)
    }
    for (m in 2:NE) {
      muT[g, m] <- muT_raw[g, m]
    }
  }

   # Sort first dimension to avoid label switching
   muT[1:G, 1] <- sort(muT_raw[1:G, 1])

   # Prior on genotype effect
  for(i in 1:I) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  # # Prior on variance matrix
  # st_inv <- inverse(st)
  # st[1,1] <- st1^2
  # st[2,2] <- st2^2
  # st[1,2] <- rho * st1 * st2
  # st[2,1] <- st[1,2]
  #
  # st1 ~ dt(0, 10^-2, 1)T(0,)
  # st2 ~ dt(0, 10^-2, 1)T(0,)
  # rho ~ dunif(-1, 1)

  st_inv ~ dwish(stE, NE)
  st <- inverse(st_inv)

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
}
"


# Set up the data
model_data <- list(N = N, y = dat$df$y, G = G, I = I, NE = NE, gen = dat$df$gen,
                   t = dat$t, stE = diag(NE))#, alpha = rep(1,G))

# Choose the parameters to watch
model_parameters <- c("g", "e", "env", "pi", "muT", "mu", "sigma", "st")
#model_parameters <- c("g")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)


# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
plot(model_run)
print(model_run)


# Plot the posterior cluster membership
qplot(model_run$BUGSoutput$median$env, dat$df$env) +
  geom_jitter(width = 0.1, height = 0.1) + theme_light() +
  labs(x = 'environment cluster', y = 'environment cluster estimated')

# Overall predictions
qplot(model_run$BUGSoutput$median$mu, dat$df$y) +
  geom_abline() + theme_light() +
  labs(x = expression(mu), y = 'y')

# Prediction of genotype effects
qplot(model_run$BUGSoutput$median$g, dat$g) +
  geom_abline() + theme_light() +
  labs(x = 'genotype estimated', y = 'genotype')

# Prediction of environment effects
qplot(model_run$BUGSoutput$median$e, dat$e) +
  geom_abline() + theme_light() +
  labs(x = 'environment estimated', y = 'environment')


## -------------------------------------------------------------------------- ##
## ------------------ Multivariate model (gen and env) --------------------- ##
## -------------------------------------------------------------------------- ##


rm(list = ls())

genData <- function(N, G1, G2, NG, NE, muG, muE, sg, se, sgroup, sigma){

  g <- rnorm(G1, 0, sg)
  e <- rnorm(G2, 0, se)

  # cluster gen
  gen <- rep(NA, N)
  g_group <- matrix(NA, ncol = NG, nrow = N)
  for (i in 1:N) {
    gen[i] <- sample(1:G1, size = 1, replace = TRUE)
    g_group[i, ] <- rmvnorm(1, muG[gen[i], ], sgroup)
  }

  # cluster env
  env <- rep(NA, N)
  e_group <- matrix(NA, ncol = NE, nrow = N)
  for (i in 1:N) {
    env[i] <- sample(1:G2, size = 1, replace = TRUE)
    e_group[i, ] <- rmvnorm(1, muE[env[i], ], sgroup)
  }

  # generate mean of the model
  mu_y <- g[gen] + e[env]
  y <- rnorm(N, mu_y, sigma)

  df <- data.frame(y = y, env = env, gen = gen)

  return(list(df = df, g = g, e = e, g_group = g_group, e_group = e_group))

}


G1 <- 3
G2 <- 2
N <- 200
NE <- 2
NG <- 2
muG <- matrix(c(-5, 0, 2, -5, 0, 2), ncol = 2, nrow = G1)
muE <- matrix(c(-5, 5, 0, 10), ncol = 2, nrow = G2)
#muE <- matrix(c(-5, 5, 0, 10), ncol = 2, nrow = G2)
sg <- 1
se <- 1
sgroup <-  diag(2)
sigma <- 1

set.seed(02)
dat <- genData(N = N, G1 = G1, G2 = G2, NG = NG, NE = NE, muG = muG, muE = muE, sg = sg, se = se, sgroup = sgroup, sigma = sigma)

plot(dat$g_group)
plot(dat$e_group)


# Jags code to fit the model to the simulated data
model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] #gen and env are parameters

    # Clustering model
    gen[i] ~ dcat(piG[i, 1:G1])
    env[i] ~ dcat(piE[i, 1:G2])

    # Continuous variables
    g_group[i, 1:NG] ~ dmnorm(mu_gen[gen[i], 1:NG], sgroup_inv_gen)
    e_group[i, 1:NE] ~ dmnorm(mu_env[env[i], 1: NE], sgroup_inv_env)

    for (g in 1:G1) {
      exp_theta_gen[i, g] <- exp(theta_gen[i, g])
      piG[i, g] <- exp(theta_gen[i, g]) / sum(exp_theta_gen[i, 1:G1])
      theta_gen[i, g] ~ dnorm(0, 6^-2)
    }


    for (g in 1:G2) {
      exp_theta_env[i, g] <- exp(theta_env[i, g])
      piE[i, g] <- exp(theta_env[i, g]) / sum(exp_theta_env[i, 1:G2])
      theta_env[i, g] ~ dnorm(0, 6^-2)
    }

   }

  # Priors on means gen
  for (g in 1:G1) {
    for (m in 1:NG) {
      mu_gen_raw[g, m] ~ dnorm(0, 100^-2)
    }
    for (m in 2:NG) {
      mu_gen[g, m] <- mu_gen_raw[g, m]
    }
  }

   # Sort first dimension to avoid label switching
   mu_gen[1:G1, 1] <- sort(mu_gen_raw[1:G1, 1])



   # Priors on means env
  for (g in 1:G2) {
    for (m in 1:NE) {
      mu_env_raw[g, m] ~ dnorm(0, 100^-2)
    }
    for (m in 2:NE) {
      mu_env[g, m] <- mu_env_raw[g, m]
    }
  }

   # Sort first dimension to avoid label switching
   mu_env[1:G2, 1] <- sort(mu_env_raw[1:G2, 1])

   # Prior on genotype effect
  for(i in 1:G1) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G2) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  # # Prior on variance matrix
  # sgroup_inv <- inverse(sgroup)
  # sgroup[1,1] <- sgroup1^2
  # sgroup[2,2] <- sgroup2^2
  # sgroup[1,2] <- rho * sgroup1 * sgroup2
  # sgroup[2,1] <- sgroup[1,2]
  #
  # sgroup1 ~ dt(0, 10^-2, 1)T(0,)
  # sgroup2 ~ dt(0, 10^-2, 1)T(0,)
  # rho ~ dunif(-1, 1)

  sgroup_inv_gen ~ dwish(sG, NG)
  sgroup_gen <- inverse(sgroup_inv_gen)

  sgroup_inv_env ~ dwish(sE, NE)
  sgroup_env <- inverse(sgroup_inv_env)

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
}
"


# Set up the data
model_data <- list(N = N, y = dat$df$y, G1 = G1, G2 = G2, NG = NG, NE = NE, gen = dat$df$gen, env = dat$df$env,
                   g_group = dat$g_group, e_group = dat$e_group, sG = diag(NG), sE = diag(NE)) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
str(model_data)
# Choose the parameters to watch
model_parameters <- c("g", "e", "gen", "env", "piG", "piE", "mu_env", "mu_gen",
                      "mu", "sigma", "sgroup_gen", 'sgroup_env')

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
plot(model_run)
print(model_run)


# Plot the posterior cluster membership
gen_c <- qplot(model_run$BUGSoutput$median$gen, dat$df$gen) +
  geom_jitter(width = 0.1, height = 0.1) + theme_light() +
  labs(x = 'genotype cluster', y = 'genotype cluster estimated')


env_c <- qplot(model_run$BUGSoutput$median$env, dat$df$env) +
  geom_jitter(width = 0.1, height = 0.1) + theme_light() +
  labs(x = 'environment cluster', y = 'environment cluster estimated')

gen_c + env_c

# Overall predictions
qplot(model_run$BUGSoutput$median$mu, dat$df$y) +
  geom_abline() + theme_light() +
  labs(x = expression(mu), y = 'y')

# Prediction of genotype effects
gen_p <- qplot(model_run$BUGSoutput$median$g, dat$g) +
  geom_abline() + theme_light() +
  labs(x = 'genotype estimated', y = 'genotype')

# Prediction of environment effects
env_p <- qplot(model_run$BUGSoutput$median$e, dat$e) +
  geom_abline() + theme_light() +
  labs(x = 'environment estimated', y = 'environment')

gen_p + env_p


