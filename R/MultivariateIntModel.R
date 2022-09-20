library(R2jags)
library(mvtnorm) # for multivariate normal distribution
library(ggplot2)

## -------------------------------------------------------------------------- ##
## ------------------ Multivariate interaction model  (env) --------------------- ##
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
