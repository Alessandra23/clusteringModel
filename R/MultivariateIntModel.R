library(R2jags)
library(mvtnorm) # for multivariate normal distribution
library(ggplot2)
library(patchwork)

## -------------------------------------------------------------------------- ##
## ------------------ Multivariate interaction model  (env) --------------------- ##
## -------------------------------------------------------------------------- ##


rm(list = ls())

generateBlin <- function(index, Q, stheta = 1){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  sqrtTheta <-  vector()
  for(q in 1:Q){
    for (i in 1:index) {
      theta[i,q] <- rnorm(1, 0, stheta)
    }
    #theta[index, q] <- -sum(theta[1:(index-1), q])
    m <- apply(theta, 2, mean)
    thetaN <- as.matrix(apply(theta, 1, function(x){x-m}))

    if(Q>1){
      thetaN <- t(thetaN)
    }

    sqrtTheta[q] <- sqrt(1/sum(thetaN[,q]^2))
    for (i in 1:index) {
      variable[i,q] <- (thetaN[i,q])*sqrtTheta[q]
    }
  }

  return(variable)
}

genData <- function(N, G, I, NE, muT, sigma, sg, se, lambda){

  Q <- length(lambda)
  g <- rnorm(I, 0, sg)
  e <- rnorm(G, 0, se)
  gen <- rep(1:I, each = N/I)

  st <- MCMCpack::riwish(NE, diag(NE))

  # cluster membership
  env <- rep(NA, N)
  t <- matrix(NA, ncol = NE, nrow = N)
  for (i in 1:N) {
    env[i] <- sample(1:G, size = 1, replace = TRUE)
    t[i, ] <- mvtnorm::rmvnorm(1, muT[env[i], ], st)
  }
  # env <- sample(1:G, size = N, replace = TRUE)
  # t <- rnorm(N, muT[env], st)

  # generate lambda, gamma, delta and kappa
  gamma <- generateBlin(I, Q)
  delta <- generateBlin(G, Q)

  # generate bilinear term
  blin <- rep(0, N)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k] * gamma[gen, k] * delta[env, k]
  }

  # generate mean of the model
  mu_y <- g[gen] + e[env] + blin
  y <- rnorm(N, mu_y, sigma)

  df <- data.frame(y = y, env = env, gen = gen, blin = blin)

  return(list(df = df, g = g, e = e, t = t, Q = Q))

}


G <- 2
N <- 90
NE <- 2
muT <- matrix(c(-5, 2, -5, 2), ncol = 2, nrow = G)
#muT <- matrix(c(-5, 0, 2, -5, 0, 2), ncol = 2, nrow = G)
#muT <- matrix(c(0.025, 0.075, 0.125, 30,42,50), ncol = 2, nrow = G)
sigma <- 1
sg <- 10
se <- 10
I <- 9
lambda <- 10

set.seed(02)
dat <- genData(N = N, G = G, NE = NE, muT = muT, sigma = sigma, I = I, sg = sg, se = se, lambda = lambda)
plot(dat$t)


# Jags code to fit the model to the simulated data
model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] + blin[i] #env here is a parameter but gen is not
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

    # # Clustering model
    t[i, 1:NE] ~ dmnorm(muT[env[i], 1:NE], st_inv)
    env[i] ~ dcat(pi[i, 1:G])

    for (j in 1:G) {
      exp_theta[i, j] <- exp(theta[i, j])
      pi[i, j] <- exp(theta[i, j]) / sum(exp_theta[i, 1:G])
      theta[i, j] ~ dnorm(0, 6^-2)
    }

   }

  # # Prior on cluster membership
  # for (g in 1:G) {
  #     pi[1:NE, g] ~ ddirch(alpha)
  #   }


  # Priors on means
  for (j in 1:G) {
    for (m in 1:NE) {
      muT_raw[j, m] ~ dnorm(0, 100^-2)
    }
    for (m in 2:NE) {
      muT[j, m] <- muT_raw[j, m]
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

  # Priors on gamma
  for(q in 1:Q){
    for(i in 1:I){
      thetaG[i,q] ~ dnorm(0,1)
    }
    mG[q] = sum(thetaG[1:I,q])/I
    for(i in 1:I){
    thetaGNew[i,q] = thetaG[i,q] - mG[q]
    }
    sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:I,q]^2 + 0.000001)))
    for(i in 1:I){
      gamma[i,q] = thetaGNew[i,q]*sqrtThetaG[q]
    }
  }

   # Priors on delta
   for(q in 1:Q){
    for(j in 1:G){
      thetaD[j,q] ~ dnorm(0,1)
    }
    mD[q] = sum(thetaD[1:G,q])/G
    for(j in 1:G){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:G,q]^2+ 0.000001)))
    for(j in 1:G){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, 100^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  st_inv ~ dwish(stE, NE)
  st <- inverse(st_inv)

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
}
"


# Set up the data
Gused <- 10
model_data <- list(N = N, y = dat$df$y, G = Gused, I = I, NE = NE, gen = dat$df$gen,
                   t = dat$t, stE = diag(NE), Q = dat$Q)#, alpha = rep(1,G))

# Choose the parameters to watch
model_parameters <- c("g", "e", "env", "pi", "muT", "mu", "sigma", "st", 'blin')
#model_parameters <- c("g")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)


plot(model_run)
print(model_run)

## see rmse
round(caret::RMSE(model_run$BUGSoutput$median$mu, dat$df$y),4)
set.seed(04)
dat2 <- genData(N = N, G = G, NE = NE, muT = muT, sigma = sigma, I = I, sg = sg, se = se, lambda = lambda)
round(caret::RMSE(model_run$BUGSoutput$median$mu, dat2$df$y),4)


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

# Prediction of bilinear effect
qplot(model_run$BUGSoutput$median$blin, dat$df$blin) +
  geom_abline() + theme_light() +
  labs(x = expression(hat(int)), y = 'int')



## --------------------------------------------------------------------------------------- ##
## ------------------ Multivariate interaction model  (gen and env) --------------------- ##
## -------------------------------------------------------------------------------------- ##


rm(list = ls())

generateBlin <- function(index, Q, stheta = 1){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  sqrtTheta <-  vector()
  for(q in 1:Q){
    for (i in 1:index) {
      theta[i,q] <- rnorm(1, 0, stheta)
    }
    #theta[index, q] <- -sum(theta[1:(index-1), q])
    m <- apply(theta, 2, mean)
    thetaN <- as.matrix(apply(theta, 1, function(x){x-m}))

    if(Q>1){
      thetaN <- t(thetaN)
    }

    sqrtTheta[q] <- sqrt(1/sum(thetaN[,q]^2))
    for (i in 1:index) {
      variable[i,q] <- (thetaN[i,q])*sqrtTheta[q]
    }
  }

  return(variable)
}

genData <- function(N, G1, G2, NG, NE, muG, muE, sg, se, sigma, lambda){

  Q <- length(lambda)
  g <- rnorm(G1, 0, sg)
  e <- rnorm(G2, 0, se)

  sgroupG <- MCMCpack::riwish(NG, diag(NG))
  sgroupE <- MCMCpack::riwish(NE, diag(NE))

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

  # generate lambda, gamma, delta and kappa
  gamma <- generateBlin(G1, Q)
  delta <- generateBlin(G2, Q)

  # generate bilinear term
  blin <- rep(0, N)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k] * gamma[gen, k] * delta[env, k]
  }

  # generate mean of the model
  mu_y <- g[gen] + e[env] + blin
  y <- rnorm(N, mu_y, sigma)

  df <- data.frame(y = y, env = env, gen = gen, blin = blin)

  return(list(df = df, g = g, e = e, g_group = g_group, e_group = e_group, Q = Q))

}



G1 <- 3
G2 <- 3
N <- 90
NE <- 2
NG <- 2
#muG <- matrix(c(-5, 0, 2, -5, 0, 2), ncol = NG, nrow = G1)
muG <- matrix(runif(NG*G1, -50,50), ncol = NG, nrow = G1)
#muG <- matrix(c(-5, 5, 0, 10), ncol = 2, nrow = G1)
#muE <- matrix(c(-5, 5, 0, 10), ncol = 2, nrow = G2)
muE <- matrix(runif(NE*G2, -50,50), ncol = NE, nrow = G2)
#muE <- matrix(c(-5, 0, 2, -5, 0, 2), ncol = NE, nrow = G2)
sg <- 1
se <- 1
#sgroup <-  diag(2)
sigma <- 1
lambda <- 100

set.seed(02)
dat <- genData(N = N, G1 = G1, G2 = G2, NG = NG, NE = NE,
               muG = muG, muE = muE, sg = sg, se = se,
               sigma = sigma, lambda = lambda)

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
    mu[i] = g[gen[i]] + e[env[i]] + blin[i] #gen and env are parameters
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

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

  # Priors on gamma
  for(q in 1:Q){
    for(i in 1:G1){
      thetaG[i,q] ~ dnorm(0,1)
    }
    mG[q] = sum(thetaG[1:G1,q])/G1
    for(i in 1:G1){
    thetaGNew[i,q] = thetaG[i,q] - mG[q]
    }
    sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:G1,q]^2 + 0.000001)))
    for(i in 1:G1){
      gamma[i,q] = thetaGNew[i,q]*sqrtThetaG[q]
    }
  }

   # Priors on delta
   for(q in 1:Q){
    for(j in 1:G2){
      thetaD[j,q] ~ dnorm(0,1)
    }
    mD[q] = sum(thetaD[1:G2,q])/G2
    for(j in 1:G2){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:G2,q]^2+ 0.000001)))
    for(j in 1:G2){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, 100^-2)T(0,)
  }
  lambda = sort(lambda_raw)

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
                   g_group = dat$g_group, e_group = dat$e_group, sG = diag(NG), sE = diag(NE), Q = dat$Q) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
str(model_data)
# Choose the parameters to watch
model_parameters <- c("g", "e", "gen", "env", "piG", "piE", "mu_env", "mu_gen",
                      "mu", "sigma", "sgroup_gen", 'sgroup_env', 'blin')

# Run the model
model_run <- R2jags::jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
plot(model_run)
print(model_run)


# Plot the posterior cluster membership
p1 <- qplot(model_run$BUGSoutput$median$env, dat$df$env) +
  geom_jitter(width = 0.1, height = 0.1) + theme_light() +
  labs(x = 'environment cluster', y = 'environment cluster estimated')

p2 <- qplot(model_run$BUGSoutput$median$gen, dat$df$gen) +
  geom_jitter(width = 0.1, height = 0.1) + theme_light() +
  labs(x = 'genotypic cluster', y = 'genotypic cluster estimated')


p1 + p2

# Overall predictions
qplot(model_run$BUGSoutput$median$mu, dat$df$y) +
  geom_abline() + theme_light() +
  labs(x = expression(hat(y)), y = 'y')
caret::RMSE(model_run$BUGSoutput$median$mu, dat$df$y)

# Prediction of genotype effects
qplot(model_run$BUGSoutput$median$g, dat$g) +
  geom_abline() + theme_light() +
  labs(x = 'genotype estimated', y = 'genotype')

# Prediction of environment effects
qplot(model_run$BUGSoutput$median$e, dat$e) +
  geom_abline() + theme_light() +
  labs(x = 'environment estimated', y = 'environment')

qplot(model_run$BUGSoutput$median$blin, dat$df$blin) +
  geom_abline() + theme_light() +
  labs(x = expression(hat(int)), y = 'int')


model_run$BUGSoutput$median$blin
dat$df$blin



# Clustering + BAMMIT -----------------------------------------------------

model_code_clustG <- "
model
{
  # Likelihood
   for (i in 1:N) {
    gen[i] ~ dcat(piG[i, 1:G1])

    # Continuous variables
    g_group[i, 1:NG] ~ dmnorm(mu_gen[gen[i], 1:NG], sgroup_inv_gen)

    for (j in 1:G1) {
      exp_theta_gen[i, j] <- exp(theta_gen[i, j])
      piG[i, j] <- exp(theta_gen[i, j]) / sum(exp_theta_gen[i, 1:G1])
      theta_gen[i, j] ~ dnorm(0, 6^-2)
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

  sgroup_inv_gen ~ dwish(sG, NG)
  sgroup_gen <- inverse(sgroup_inv_gen)

}
"

model_code_clustE <- "
model
{
  # Likelihood
  for (i in 1:N) {
    env[i] ~ dcat(piE[i, 1:G2])

    # Continuous variables
    e_group[i, 1:NE] ~ dmnorm(mu_env[env[i], 1: NE], sgroup_inv_env)

    for (j in 1:G2) {
      exp_theta_env[i, j] <- exp(theta_env[i, j])
      piE[i, j] <- exp(theta_env[i, j]) / sum(exp_theta_env[i, 1:G2])
      theta_env[i, j] ~ dnorm(0, 6^-2)
    }

   }

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

  sgroup_inv_env ~ dwish(sE, NE)
  sgroup_env <- inverse(sgroup_inv_env)
}
"



# Set up the data
usedG <- 2
usedE <- 2
model_data_groupsG <- list(N = N, G1 = G1, NG = NG, gen = dat$df$gen,
                          g_group = dat$g_group, sG = diag(NG))
model_data_groupsE <- list(N = N, G2 = G2, NE = NE, env = dat$df$env,
                           e_group = dat$e_group, sE = diag(NE))

# Choose the parameters to watch
model_parameters_groupsG <- c("gen", "piG","mu_gen", "sgroup_gen")
model_parameters_groupsE <- c("env", "piE", "mu_env", 'sgroup_env')

# Run the model
model_run_groupsG <- jags(
  data = model_data_groupsG,
  parameters.to.save = model_parameters_groupsG,
  model.file = textConnection(model_code_clustG)
)

model_run_groupsE <- jags(
  data = model_data_groupsE,
  parameters.to.save = model_parameters_groupsE,
  model.file = textConnection(model_code_clustE)
)


groupsG <- model_run_groupsG$BUGSoutput$median$gen
groupsE <- model_run_groupsE$BUGSoutput$median$env
as.integer(groupsG)
as.integer(groupsE)


model_code_BAMMIT <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] + blin[i] #gen and env are parameters
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

   }


   # Prior on genotype effect
  for(i in 1:G1) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G2) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  # Priors on gamma
  for(q in 1:Q){
    for(i in 1:G1){
      thetaG[i,q] ~ dnorm(0,1)
    }
    mG[q] = sum(thetaG[1:G1,q])/G1
    for(i in 1:G1){
    thetaGNew[i,q] = thetaG[i,q] - mG[q]
    }
    sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:G1,q]^2 + 0.000001)))
    for(i in 1:G1){
      gamma[i,q] = thetaGNew[i,q]*sqrtThetaG[q]
    }
  }

   # Priors on delta
   for(q in 1:Q){
    for(j in 1:G2){
      thetaD[j,q] ~ dnorm(0,1)
    }
    mD[q] = sum(thetaD[1:G2,q])/G2
    for(j in 1:G2){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:G2,q]^2+ 0.000001)))
    for(j in 1:G2){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, 100^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
}
"


# Set up the data
model_data_BAMMIT <- list(N = N, y = dat$df$y, G1 = G1, G2 = G2, gen = groupsG, env = groupsE, Q = dat$Q) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
str(model_data_BAMMIT)
# Choose the parameters to watch
model_parameters_BAMMIT <- c("g", "e", "mu", "sigma",  'blin')

# Run the model
model_run_BAMMIT <- R2jags::jags(
  data = model_data_BAMMIT,
  parameters.to.save = model_parameters_BAMMIT,
  model.file = textConnection(model_code_BAMMIT)
)

plot(model_run)
print(model_run)

# Overall predictions
caret::RMSE(model_run$BUGSoutput$median$mu, dat$df$y)














# Compare both methods ----------------------------------------------------

rm(list = ls())


generateBlin <- function(index, Q, stheta = 1){
  theta <- matrix(NA, nrow = index, ncol = Q)
  variable <- matrix(NA, nrow = index, ncol = Q)
  sqrtTheta <-  vector()
  for(q in 1:Q){
    for (i in 1:index) {
      theta[i,q] <- rnorm(1, 0, stheta)
    }
    #theta[index, q] <- -sum(theta[1:(index-1), q])
    m <- apply(theta, 2, mean)
    thetaN <- as.matrix(apply(theta, 1, function(x){x-m}))

    if(Q>1){
      thetaN <- t(thetaN)
    }

    sqrtTheta[q] <- sqrt(1/sum(thetaN[,q]^2))
    for (i in 1:index) {
      variable[i,q] <- (thetaN[i,q])*sqrtTheta[q]
    }
  }

  return(variable)
}

genData <- function(N, G1, G2, NG, NE, muG, muE, sg, se, sigma, lambda){

  Q <- length(lambda)
  g <- rnorm(G1, 0, sg)
  e <- rnorm(G2, 0, se)

  sgroupG <- MCMCpack::riwish(NG, diag(NG))
  sgroupE <- MCMCpack::riwish(NE, diag(NE))

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

  # generate lambda, gamma, delta and kappa
  gamma <- generateBlin(G1, Q)
  delta <- generateBlin(G2, Q)

  # generate bilinear term
  blin <- rep(0, N)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k] * gamma[gen, k] * delta[env, k]
  }

  # generate mean of the model
  mu_y <- g[gen] + e[env] + blin
  y <- rnorm(N, mu_y, sigma)

  df <- data.frame(y = y, env = env, gen = gen, blin = blin)

  return(list(df = df, g = g, e = e, g_group = g_group, e_group = e_group, Q = Q))

}



simData <- function(N = 120, trueC, NE = 3, NG = 3, seedTrain = 02, seedTest = 04){

  G1 <- trueC
  G2 <- trueC
  N <- N
  NE <- NE
  NG <- NG
  sg <- 1
  se <- 1
  #sgroup <-  diag(2)
  sigma <- 1
  lambda <- 100

  set.seed(seedTrain)
  muG <- matrix(runif(NG*G1, -50,50), ncol = NG, nrow = G1)
  muE <- matrix(runif(NE*G2, -50,50), ncol = NE, nrow = G2)
  datTrain <- genData(N = N, G1 = G1, G2 = G2, NG = NG, NE = NE,
                 muG = muG, muE = muE, sg = sg, se = se,
                 sigma = sigma, lambda = lambda)

  set.seed(seedTest)
  muG <- matrix(runif(NG*G1, -50,50), ncol = NG, nrow = G1)
  muE <- matrix(runif(NE*G2, -50,50), ncol = NE, nrow = G2)
  datTest <- genData(N = N, G1 = G1, G2 = G2, NG = NG, NE = NE,
                      muG = muG, muE = muE, sg = sg, se = se,
                      sigma = sigma, lambda = lambda)

  return(list(datTrain = datTrain,
              datTest = datTest))

}


models <- function(dat, usedC, N = 120, NE = 3, NG = 3){

  # CBAMMIT

  model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] + blin[i] #gen and env are parameters
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

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

  # Priors on gamma
  for(q in 1:Q){
    for(i in 1:G1){
      thetaG[i,q] ~ dnorm(0,1)
    }
    mG[q] = sum(thetaG[1:G1,q])/G1
    for(i in 1:G1){
    thetaGNew[i,q] = thetaG[i,q] - mG[q]
    }
    sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:G1,q]^2 + 0.000001)))
    for(i in 1:G1){
      gamma[i,q] = thetaGNew[i,q]*sqrtThetaG[q]
    }
  }

   # Priors on delta
   for(q in 1:Q){
    for(j in 1:G2){
      thetaD[j,q] ~ dnorm(0,1)
    }
    mD[q] = sum(thetaD[1:G2,q])/G2
    for(j in 1:G2){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:G2,q]^2+ 0.000001)))
    for(j in 1:G2){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, 100^-2)T(0,)
  }
  lambda = sort(lambda_raw)

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
  model_data <- list(N = N, y = dat$df$y, G1 = usedC, G2 = usedC, NG = NG, NE = NE, gen = dat$df$gen, env = dat$df$env,
                     g_group = dat$g_group, e_group = dat$e_group, sG = diag(NG), sE = diag(NE), Q = dat$Q) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
  str(model_data)
  # Choose the parameters to watch
  model_parameters <- c("g", "e", "gen", "env", "piG", "piE", "mu_env", "mu_gen",
                        "mu", "sigma", "sgroup_gen", 'sgroup_env', 'blin')

  # Run the model
  model_run <- R2jags::jags(
    data = model_data,
    parameters.to.save = model_parameters,
    model.file = textConnection(model_code)
  )


  # C + BAMMIT

  model_code_clustG <- "
model
{
  # Likelihood
   for (i in 1:N) {
    gen[i] ~ dcat(piG[i, 1:G1])

    # Continuous variables
    g_group[i, 1:NG] ~ dmnorm(mu_gen[gen[i], 1:NG], sgroup_inv_gen)

    for (j in 1:G1) {
      exp_theta_gen[i, j] <- exp(theta_gen[i, j])
      piG[i, j] <- exp(theta_gen[i, j]) / sum(exp_theta_gen[i, 1:G1])
      theta_gen[i, j] ~ dnorm(0, 6^-2)
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

  sgroup_inv_gen ~ dwish(sG, NG)
  sgroup_gen <- inverse(sgroup_inv_gen)

}
"

  model_code_clustE <- "
model
{
  # Likelihood
  for (i in 1:N) {
    env[i] ~ dcat(piE[i, 1:G2])

    # Continuous variables
    e_group[i, 1:NE] ~ dmnorm(mu_env[env[i], 1: NE], sgroup_inv_env)

    for (j in 1:G2) {
      exp_theta_env[i, j] <- exp(theta_env[i, j])
      piE[i, j] <- exp(theta_env[i, j]) / sum(exp_theta_env[i, 1:G2])
      theta_env[i, j] ~ dnorm(0, 6^-2)
    }

   }

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

  sgroup_inv_env ~ dwish(sE, NE)
  sgroup_env <- inverse(sgroup_inv_env)
}
"

  model_data_groupsG <- list(N = N, G1 = usedC, NG = NG, gen = dat$df$gen,
                             g_group = dat$g_group, sG = diag(NG))
  model_data_groupsE <- list(N = N, G2 = usedC, NE = NE, env = dat$df$env,
                             e_group = dat$e_group, sE = diag(NE))

  # Choose the parameters to watch
  model_parameters_groupsG <- c("gen", "piG","mu_gen", "sgroup_gen")
  model_parameters_groupsE <- c("env", "piE", "mu_env", 'sgroup_env')

  # Run the model
  model_run_groupsG <- jags(
    data = model_data_groupsG,
    parameters.to.save = model_parameters_groupsG,
    model.file = textConnection(model_code_clustG)
  )

  model_run_groupsE <- jags(
    data = model_data_groupsE,
    parameters.to.save = model_parameters_groupsE,
    model.file = textConnection(model_code_clustE)
  )


  groupsG <- as.integer(model_run_groupsG$BUGSoutput$median$gen)
  groupsE <- as.integer(model_run_groupsE$BUGSoutput$median$env)


  model_code_BAMMIT <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] + blin[i] #gen and env are parameters
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

   }


   # Prior on genotype effect
  for(i in 1:G1) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G2) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  # Priors on gamma
  for(q in 1:Q){
    for(i in 1:G1){
      thetaG[i,q] ~ dnorm(0,1)
    }
    mG[q] = sum(thetaG[1:G1,q])/G1
    for(i in 1:G1){
    thetaGNew[i,q] = thetaG[i,q] - mG[q]
    }
    sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:G1,q]^2 + 0.000001)))
    for(i in 1:G1){
      gamma[i,q] = thetaGNew[i,q]*sqrtThetaG[q]
    }
  }

   # Priors on delta
   for(q in 1:Q){
    for(j in 1:G2){
      thetaD[j,q] ~ dnorm(0,1)
    }
    mD[q] = sum(thetaD[1:G2,q])/G2
    for(j in 1:G2){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:G2,q]^2+ 0.000001)))
    for(j in 1:G2){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, 100^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
}
"


  # Set up the data
  model_data_BAMMIT <- list(N = N, y = dat$df$y, G1 = usedC, G2 = usedC, gen = groupsG, env = groupsE, Q = dat$Q) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
  #str(model_data_BAMMIT)
  # Choose the parameters to watch
  model_parameters_BAMMIT <- c("g", "e", "mu", "sigma",  'blin')

  # Run the model
  model_run_BAMMIT <- R2jags::jags(
    data = model_data_BAMMIT,
    parameters.to.save = model_parameters_BAMMIT,
    model.file = textConnection(model_code_BAMMIT)
  )


  return(list(CBAMMIT = model_run,
              ClustBammit = model_run_BAMMIT))

}


# gen data
# trueC <- c(2,10,20)
# usedC <- c(2,10,20,50)

dat <- simData(N = 120, trueC = 2, NE = 3, NG = 3, seedTrain = 02, seedTest = 04)

# ran the models
modelsResults <- vector('list')

modelsResults$usedC2 <- models(dat = dat$datTrain, usedC = 2, N = 120, NE = 3, NG = 3)
modelsResults$usedC10 <- models(dat = dat$datTrain, usedC = 10, N = 120, NE = 3, NG = 3)
modelsResults$usedC20 <- models(dat = dat$datTrain, usedC = 20, N = 120, NE = 3, NG = 3)
modelsResults$usedC50 <- models(dat = dat$datTrain, usedC = 50, N = 120, NE = 3, NG = 3)


# Train
round(caret::RMSE(modelsResults$usedC2$CBAMMIT$BUGSoutput$median$mu, dat$datTrain$df$y), 4)
round(caret::RMSE(modelsResults$usedC2$ClustBammit$BUGSoutput$median$mu, dat$datTrain$df$y), 4)

round(caret::RMSE(modelsResults$usedC10$CBAMMIT$BUGSoutput$median$mu, dat$datTrain$df$y), 4)
round(caret::RMSE(modelsResults$usedC10$ClustBammit$BUGSoutput$median$mu, dat$datTrain$df$y), 4)

round(caret::RMSE(modelsResults$usedC20$CBAMMIT$BUGSoutput$median$mu, dat$datTrain$df$y), 4)
round(caret::RMSE(modelsResults$usedC20$ClustBammit$BUGSoutput$median$mu, dat$datTrain$df$y), 4)

round(caret::RMSE(modelsResults$usedC50$CBAMMIT$BUGSoutput$median$mu, dat$datTrain$df$y), 4)
round(caret::RMSE(modelsResults$usedC50$ClustBammit$BUGSoutput$median$mu, dat$datTrain$df$y), 4)


# Test
round(caret::RMSE(modelsResults$usedC2$CBAMMIT$BUGSoutput$median$mu, dat$datTest$df$y), 4)
round(caret::RMSE(modelsResults$usedC2$ClustBammit$BUGSoutput$median$mu, dat$datTest$df$y), 4)

round(caret::RMSE(modelsResults$usedC10$CBAMMIT$BUGSoutput$median$mu, dat$datTest$df$y), 4)
round(caret::RMSE(modelsResults$usedC10$ClustBammit$BUGSoutput$median$mu, dat$datTest$df$y), 4)

round(caret::RMSE(modelsResults$usedC20$CBAMMIT$BUGSoutput$median$mu, dat$datTest$df$y), 4)
round(caret::RMSE(modelsResults$usedC20$ClustBammit$BUGSoutput$median$mu, dat$datTest$df$y), 4)

round(caret::RMSE(modelsResults$usedC50$CBAMMIT$BUGSoutput$median$mu, dat$datTest$df$y), 4)
round(caret::RMSE(modelsResults$usedC50$ClustBammit$BUGSoutput$median$mu, dat$datTest$df$y), 4)


plot(modelsResults$usedC2$ClustBammit$BUGSoutput$median$mu, dat$datTrain$df$y)




