library(R2jags)
library(ggplot2)

rm(list = ls())


# Generate data -----------------------------------------------------------

genData <- function(N, G1, G2, NG, NE, muG, muE, matrixG, matrixE, sg, se, sigma){

  g <- rnorm(G1, 0, sg)
  e <- rnorm(G2, 0, se)

  # g <- g - mean(g)
  # e <- e - mean(e)

  sgroupG <- MCMCpack::riwish(NG, matrixG)
  sgroupE <- MCMCpack::riwish(NE, matrixE)

  # sgroupG <- diag(NG)
  # sgroupE <- diag(NE)

  # cluster gen
  gen <- rep(NA, N)
  g_group <- matrix(NA, ncol = NG, nrow = N)
  for (i in 1:N) {
    gen[i] <- sample(1:G1, size = 1, replace = TRUE)
    g_group[i, ] <- mvtnorm::rmvnorm(1, muG[gen[i], ], sgroupG)
  }

  # cluster env
  env <- rep(NA, N)
  e_group <- matrix(NA, ncol = NE, nrow = N)
  for (i in 1:N) {
    env[i] <- sample(1:G2, size = 1, replace = TRUE)
    e_group[i, ] <- mvtnorm::rmvnorm(1, muE[env[i], ], sgroupE)
  }


  # generate mean of the model
  mu_y <- g[gen] + e[env]
  y <- rnorm(N, mu_y, sigma)

  df <- data.frame(y = y, env = env, gen = gen)

  return(list(df = df, g = g, e = e, g_group = g_group, e_group = e_group))

}



N <- 100
# number of clusters
G1 <- 3
G2 <- 3
# number of numerical variables
NE <- 2
NG <- 2
# means of teh clusters
muG <- matrix(runif(NG*G1, -50,50), ncol = NG, nrow = G1)
muE <- matrix(runif(NE*G2, -50,50), ncol = NE, nrow = G2)
# Scale matrix for riwish
matrixG <- diag(NG)
matrixE <- diag(NE)
# matrixG[1,2] <- 0.1
# matrixG[2,1] <- matrixG[1,2]
# matrixE[1,2] <- 0.1
# matrixE[2,1] <- matrixG[1,2]
# sd's
sg <- 1
se <- 1
sigma <- 1

set.seed(02)
dat <- genData(N = N, G1 = G1, G2 = G2, NG = NG, NE = NE,
               muG = muG, muE = muE, matrixG = matrixG,
               matrixE = matrixE, sg = sg, se = se, sigma = sigma)

plot(dat$g_group)
plot(dat$e_group)

qplot(dat$e_group[,1], dat$e_group[,2])

qplot(dat$e_group[,1], dat$e_group[,2]) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_bw(base_size = 16) +
  labs(x = 'temperature', y = 'rainfall')


# Jags code to fit the model to the simulated data
model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]]  # gen and env are parameters

    # Clustering model
    gen[i] ~ dcat(piG[i, 1:G1])
    env[i] ~ dcat(piE[i, 1:G2])

    # Continuous variables
    g_group[i, 1:NG] ~ dmnorm(mu_gen[gen[i], 1:NG], sgroup_inv_gen)
    e_group[i, 1:NE] ~ dmnorm(mu_env[env[i], 1: NE], sgroup_inv_env)

    for (j in 1:G1) {
      exp_theta_gen[i, j] <- exp(theta_gen[i, j])
      piG[i, j] <- exp(theta_gen[i, j]) / sum(exp_theta_gen[i, 1:G1])
      theta_gen[i, j] ~ dnorm(0, 6^-2)
    }

    for (j in 1:G2) {
      exp_theta_env[i, j] <- exp(theta_env[i, j])
      piE[i, j] <- exp(theta_env[i, j]) / sum(exp_theta_env[i, 1:G2])
      theta_env[i, j] ~ dnorm(0, 6^-2)
    }
   }

  # Priors on means gen
  for (j in 1:G1) {
    for (m in 1:NG) {
      mu_gen_raw[j, m] ~ dnorm(0, 100^-2)
    }

    for (m in 2:NG) {
      mu_gen[j, m] <- mu_gen_raw[j, m]
    }

  }

   # Sort first dimension to avoid label switching
   mu_gen[1:G1, 1] <- sort(mu_gen_raw[1:G1, 1])

   # Priors on means env
  for (j in 1:G2) {
    for (m in 1:NE) {
      mu_env_raw[j, m] ~ dnorm(0, 100^-2)
    }
    for (m in 2:NE) {
      mu_env[j, m] <- mu_env_raw[j, m]
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
model_data <- list(N = N, y = dat$df$y, G1 = G1, G2 = G2, NG = NG, NE = NE,
                   g_group = dat$g_group, e_group = dat$e_group,
                   sG = diag(NG), sE = diag(NE))
str(model_data)

# Parameters to watch
model_parameters <- c("g", "e", "gen", "env", "piG", "piE", "mu_env", "mu_gen",
                      "mu", "sigma", "sgroup_gen", 'sgroup_env')

# Run the model
model_run <- R2jags::jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# see convergence
plot(model_run)

R2jags::traceplot(model_run, varname = 'gen')

posterior_samples_sorted <- apply(model_run$BUGSoutput$sims.list$gen, 2, sort)

label.switching::label.switching(model_run$BUGSoutput)



model_run$BUGSoutput |> names()

extract(model_run)

# Plot the posterior cluster membership
qplot(dat$df$env, model_run$BUGSoutput$median$env) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_bw(base_size = 16) +
  labs(x = 'True cluster indicator', y = 'Estimated cluster indicator')

# Overall predictions
qplot(model_run$BUGSoutput$median$mu, dat$df$y) + geom_abline()

# Prediction of genotype effects
qplot(model_run$BUGSoutput$median$g, dat$g) + geom_abline()

# Prediction of environment effects
qplot(model_run$BUGSoutput$median$e, dat$e) + geom_abline()

