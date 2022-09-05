library(R2jags)

genData <- function(N, G, I, muT, sigma, sg){

  gen <- rnorm(I, 0, sg)
  g <- rep(1:I, each = N/I)

  # clustering Temp
  theta <- matrix(rnorm(N * G, 0, 3), ncol = G, nrow = N)
  pi <- exp(theta) / apply(exp(theta), 1, sum)
  Z <- rep(NA, N)
  for (i in 1:N) Z[i] <- sample(1:G, size = 1, prob = pi[i, ])
  mu_e <- rnorm(N, muT[Z], sigma)

  # generate mean of the model
  mu_y <- gen[g] + mu_e[Z]
  y <- rnorm(N, mu_y, sigma)
  df <- data.frame(y = y, temp = Z, gen = g, env = mu_e)

  return(list(y = y ,
              df = df,
              muT = muT,
              mu_e = mu_e))
  #blin = blin,
  #Q = Q))

}


G <- 3
N <- 60
muT <- c(-5, 10, 30)
sigma <- 1
sg <-  1
I <-  6
lambda <- 12
set.seed(01)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg)
hist(dat$mu_e, breaks = 30, freq = FALSE)
for (g in 1:G) curve(dnorm(x, mean = dat$muT[g])/G, col = g, add = TRUE)

ggplot(dat$df, aes(x = env)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  stat_function(fun = function(x) dnorm(x, mean = dat$muT[1])/G, color = "steelblue",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = dat$muT[2])/G, color = "firebrick",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = dat$muT[3])/G, color = "chartreuse4",size = 0.7) +
  theme_bw() +
  xlim(-15,35)

ggplot(dat$df, aes(x = y)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  theme_bw()


# Clustering function -----------------------------------------------------


model_code <- '
  model
  {
  # Likelihood
   for (i in 1:N) {
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu_e[i] ~ dnorm(muT[Z[i]], sigma_e^-2)
    mu[i] = gen[genotype[i]] + mu_e[Z[i]]
    Z[i] ~ dcat(pi[i, 1:G])
    for (g in 1:G) {
      exp_theta[i, g] <- exp(theta[i, g])
      pi[i, g] <- exp(theta[i, g]) / sum(exp_theta[i, 1:G])
      theta[i, g] ~ dnorm(0, 6^-2)
    }
   }

  # Priors

  for (g in 1:G) {
    muT_raw[g] ~ dnorm(0, 100^-2)
  }

  # Make sure these are in order to avoid label switching
  muT <- sort(muT_raw[1:G])


  # Prior on genotype effect
  for(i in 1:I) {
  gen[i] ~ dnorm(0, 100^-2) # Prior on genotype effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$y, G = G, I = I, genotype = dat$df$gen)

# Choose the parameters to watch
model_parameters <- c("muT", "sigma", "Z", "pi", "mu")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)


# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
plot(model_run)
print(model_run)




# Model 2 -----------------------------------------------------------------

library(R2jags)

genData <- function(N, G1, G2, muG, muT, sigma){

  # clustering Gen
  theta_gen <- matrix(rnorm(N * G1, 0, 1), ncol = G1, nrow = N)
  pi_gen <- exp(theta_gen) / apply(exp(theta_gen), 1, sum)
  Z_gen <- rep(NA, N)
  for (i in 1:N) Z_gen[i] <- sample(1:G1, size = 1, prob = pi_gen[i, ])
  mu_gen <- rnorm(N, muG[Z_gen], sigma_gen)

  # clustering Temp
  theta <- matrix(rnorm(N * G2, 0, 3), ncol = G2, nrow = N)
  pi <- exp(theta) / apply(exp(theta), 1, sum)
  Z_env <- rep(NA, N)
  for (i in 1:N) Z_env[i] <- sample(1:G2, size = 1, prob = pi[i, ])
  mu_e <- rnorm(N, muT[Z_env], sigma_env)

  # generate mean of the model
  mu_y <- mu_gen[Z_gen] + mu_e[Z_env]
  y <- rnorm(N, mu_y, sigma)
  df <- data.frame(y = y, Z_gen = Z_gen, Z_env = Z_env, mu_gen = mu_gen, mu_e = mu_e)

  return(list(y = y ,
              df = df))
  #blin = blin,
  #Q = Q))

}


G1 <- 3
G2 <- 3
N <- 100
muT <- c(-5, 10, 30)
muG <- c(-1, 0, 1)
sigma_gen <- 1
sigma_env <- 1
sg <-  1
I <-  6
lambda <- 12
set.seed(01)
dat <- genData(N = N, G1 = G1, G2 = G2, muG = muG, muT = muT, sigma = sigma)
hist(dat$df$mu_e, breaks = 30, freq = FALSE)
for (g in 1:G2) curve(dnorm(x, mean = muT[g])/G2, col = g, add = TRUE)

hist(dat$df$mu_gen, breaks = 10, freq = FALSE)
for (g in 1:G1) curve(dnorm(x, mean = muG[g])/G1, col = g, add = TRUE)

hist(dat$y)
library(ggplot2)

ggplot(dat$df, aes(x = mu_e)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  stat_function(fun = function(x) dnorm(x, mean = muT[1])/G, color = "steelblue",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muT[2])/G, color = "firebrick",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muT[3])/G, color = "chartreuse4",size = 0.7) +
  theme_bw() +
  xlim(-15,35)




ggplot(dat$df, aes(x = mu_gen)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  stat_function(fun = function(x) dnorm(x, mean = muG[1])/G, color = "steelblue",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muG[2])/G, color = "firebrick",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muG[3])/G, color = "chartreuse4",size = 0.7) +
  theme_bw()



ggplot(dat$df, aes(x = y)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  theme_bw()

# Clustering function -----------------------------------------------------


model_code <- '
  model
  {
  # Likelihood
   for (i in 1:N) {
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu_e[i] ~ dnorm(muT[Z_gen[i]], sigma_e^-2)
    mu_g[i] ~ dnorm(muG[Z_env[i]], sigma_g^-2)
    mu[i] = mu_g[Z_gen[i]] + mu_e[Z_env[i]]
    Z_gen[i] ~ dcat(pi_gen[i, 1:G1])
    Z_env[i] ~ dcat(pi_env[i, 1:G2])
    for (g in 1:G1) {
      exp_theta_gen[i, g] <- exp(theta_gen[i, g])
      pi_gen[i, g] <- exp(theta_gen[i, g]) / sum(exp_theta_gen[i, 1:G1])
      theta_gen[i, g] ~ dnorm(0, 6^-2)
    }

    for (g in 1:G2) {
      exp_theta_env[i, g] <- exp(theta_env[i, g])
      pi_env[i, g] <- exp(theta_env[i, g]) / sum(exp_theta_env[i, 1:G2])
      theta_env[i, g] ~ dnorm(0, 6^-2)
    }


   }

  # Priors

  for (g in 1:G1) {
    muG_raw[g] ~ dnorm(0, 100^-2)
  }

   for (g in 1:G2) {
    muT_raw[g] ~ dnorm(0, 100^-2)
  }

  # Make sure these are in order to avoid label switching
  muG <- sort(muT_raw[1:G1])
  muT <- sort(muT_raw[1:G2])

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$y, G1 = G1, G2 = G2)

# Choose the parameters to watch
model_parameters <- c("muT", "muG", "sigma", "Z_gen", "Z_env", "pi_gen", "pi_env","mu")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)


# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
plot(model_run)
print(model_run)





