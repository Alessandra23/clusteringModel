library(R2jags)
library(mvtnorm) # for multivariate normal distribution
library(ggplot2)
library(patchwork)

rm(list = ls())

# Function to meet the constraints on the interaction term

generateInt <- function(index = 10, Q = 1, stheta = 1) {

  # Generate theta matrix
  theta <- matrix(rnorm(index * Q, 0, stheta), nrow = index, ncol = Q)

  # Calculate the column means of theta
  m <- colMeans(theta)

  # Center the columns of theta
  thetaN <- sweep(theta, 2, m)

  # Calculate sqrtTheta
  sqrtTheta <- sapply(1:Q, function(q) {
    sqrt(1 / sum(thetaN[, q]^2))
  })

  # Calculate variable
  variable <- sweep(thetaN, 2, sqrtTheta, "*")

  return(variable)
}

# Model 1

genData <- function(N = 50, G1 = 4, G2 = 3, NE = 1, muE = c(-2, 5), sg = 1, se = 1, sigma = 1, seed = 02){

  set.seed(seed)

  g <- rnorm(G1, 0, sg)
  g <- g - mean(g)
  e <- rnorm(G2, 0, se)
  e <- e - mean(e)

  gen <- sample(1:G1, N, replace = TRUE)

  # cluster env
  muE <- matrix(muE, ncol = NE, nrow = G2)
  sgroupE <- MCMCpack::riwish(NE, diag(NE)) # cov matrix
  env <- rep(NA, N)
  e_group <- matrix(NA, ncol = NE, nrow = N)
  for (i in 1:N) {
    env[i] <- sample(1:G2, size = 1, replace = TRUE)
    e_group[i, ] <- mvtnorm::rmvnorm(1, muE[env[i], ], sgroupE)
  }

  # generate mean of the model
  mu_y <- g[gen] + e[env]
  y <- rnorm(N, mu_y, sigma)

  df <- data.frame(y = y, env = env, gen = gen,  g = g[gen], e = e[env], e_group = e_group)

  return(df)

}

N <- 50
NE <- 1
G1 <- 6
G2 <- 3
sg <- 1
se <- 1
sigma <- 1
muE = c(-30, 0,  30)
seed = 02

dat <- genData(N = N, G1 = G1, G2 = G2, NE = NE, muE = muE,
                sg = sg, se = se, sigma = sigma, seed = seed)

sum(unique(dat_test$e))
dat_test <- genData(N = N, G1 = G1, G2 = G2, NE = NE, muE = muE,
               sg = sg, se = se, sigma = sigma, seed = 04)

plot(dat$e_group)


# Run Jags

model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]]

    # Clustering model
    env[i] ~ dcat(piE[i,1:G2])

    # Continuous variables
    e_group[i, 1:NE] ~ dmnorm(mu_env[env[i], 1: NE], sgroup_env)

    for (j in 1:G2) {
      exp_theta_env[i, j] <- exp(theta_env[i, j])
      piE[i, j] <- exp(theta_env[i, j]) / sum(exp_theta_env[i, 1:G2])
      theta_env[i, j] ~ dnorm(0, 6^-2)
    }

   }

   # # Prior on cluster membership
   #  piE ~ ddirch(alphaE)


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
    g_initial[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G2) {
    e_initial[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  g_mean <- sum(g_initial)/G1
  e_mean <- sum(e_initial)/G2

  for(i in 1:G1) {
    g[i] <-  g_initial[i] - g_mean
  }

  for(i in 1:G2) {
    e[i] <-  e_initial[i] - e_mean
  }

  sgroup_env ~ dt(0, 10^-2, 1)T(0,)
  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
}
"


# Set up the data
usedC = 10
model_data <- list(N = N, y = dat$y, G1 = G1, G2 = usedC, NE = NE,
                   gen = dat$gen) #, alphaE = rep(1,usedC))
str(model_data)
# Choose the parameters to watch
model_parameters <- c("g", "e", "gen", "env", "piE", "mu_env", "mu", "sigma", 'sgroup_env')

# Run the model
model_run <- R2jags::jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

plot(model_run)

caret::RMSE(dat$y, model_run$BUGSoutput$mean$mu)

model = model_run$BUGSoutput

b1 <- model$mean$g
b2 <- model$mean$e
groupsG <-  dat$gen
groupsG <-  dat_test$gen
groupsE <- as.integer(model$mean$env)
#groupsE <- dat$env
yhat <- b1[groupsG] + b2[groupsE]
yhat
caret::RMSE(dat$y, yhat)


