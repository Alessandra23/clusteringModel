library(R2jags)
library(ggplot2)
library(rstan)

rm(list = ls())

# Model
# y_i = mu + b^1_{g^1_i} + b^2_{g^2_i} + ... + b^V_{g^V_i}

# V: number of variables
# g^v_i: latent variable that identifies which cluster

# Generate data -----------------------------------------------------------

genData <- function(N, B1, B2, mu_g2, sigma, sigma_b1, sigma_b2, sigma_var2_1){

  b1 <- rnorm(B1, 0, sigma_b1)
  b2 <- rnorm(B2, 0, sigma_b2)
  g1 <- rep(1:B1, each = N/B1)

  # cluster membership
  g2 <- sample(1:B2, size = N, replace = TRUE)

  # generate mean of the model
  mu_y <- b1[g1] + b2[g2]
  y <- rnorm(N, mu_y, sigma)
  var2_1 <- rnorm(N, mu_g2[g2], sigma_var2_1)
  df <- data.frame(y = y, g2 = g2, g1 = g1, var2_1 = var2_1)

  return(list(df = df, b1 = b1, b2 = b2))

}


N <- 400
B1 <- 10
B2 <- 3
mu_g2 <- c(-10,0,10)
sigma <- 1
sigma_b1 <- 1
sigma_b2 <- 1
sigma_var2_1 <- 1

set.seed(02)
dat <- genData(N = N, B1 = B1, B2 = B2, mu_g2 = mu_g2,
               sigma = sigma, sigma_b1 = sigma_b1,
               sigma_b2 = sigma_b2, sigma_var2_1 = sigma_var2_1)

# indexes for train and test
train <- sample(1:N, size = 200)
test <- (1:N)[-train]

# Run model
usedB2 <- 3
#model_code <- stan_model("~/Documents/GitHub/clusteringModel/R/model_ge.stan")
model_data <- list(N = length(train),
                   y = dat$df$y[train],
                   y_test = dat$df$y[test],
                   B2 = usedB2,
                   B1 = B1,
                   g1 = dat$df$g1[train],
                   var2_1 = dat$df$var2_1[train],
                   alpha = rep(1,usedB2),
                   N_test = length(test),
                   var2_1_test = dat$df$var2_1[test],
                   g1_test = dat$df$g1[test])
model_run <- stan(file = "~/Documents/GitHub/clusteringModel/R/model_ge.stan",
                  data = model_data,
                  chains = 4,
                  warmup = 1000,
                  iter = 2000,
                  cores = 1)
