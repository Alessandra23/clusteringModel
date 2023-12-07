library(R2jags)
library(ggplot2)
library(patchwork)

rm(list = ls())

# Model 1: b^1 + b^2 ------------------------------------------------------
# In this model the variable b^1 is categorical

# Function to generate data
genData <- function(N, B1, B2, mu_var2_1, sigma, sigma_b1, sigma_b2, sigma_var2_1) {
  b1 <- rnorm(B1, 0, sigma_b1)
  b2 <- rnorm(B2, 0, sigma_b2)
  g1 <- rep(1:B1, each = N / B1)

  # cluster membership
  g2 <- sample(1:B2, size = N, replace = TRUE)

  # generate mean of the model
  mu_y <- b1[g1] + b2[g2]
  y <- rnorm(N, mu_y, sigma)
  var2_1 <- rnorm(N, mu_var2_1[g2], sigma_var2_1)
  df <- data.frame(y = y, g2 = g2, g1 = g1, var2_1 = var2_1)

  return(list(df = df, b1 = b1, b2 = b2))
}

# set up the values
N <- 400
B1 <- 10
B2 <- 3
mu_var2_1 <- c(-10, 0, 10)
sigma <- 1
sigma_b1 <- 1
sigma_b2 <- 1
sigma_var2_1 <- 1

# simulate a data set
set.seed(02)
dat <- genData(
  N = N, B1 = B1, B2 = B2, mu_var2_1 = mu_var2_1,
  sigma = sigma, sigma_b1 = sigma_b1,
  sigma_b2 = sigma_b2, sigma_var2_1 = sigma_var2_1
)

# plot data

hist(dat$df$var2_1, breaks = 30, freq = FALSE)
for (g in 1:B2) curve(dnorm(x, mean = mu_var2_1[g], sd = sigma_var2_1) / B2, col = g, add = TRUE)

# indexes for train and test
train <- sample(1:N, size = 200)
test <- (1:N)[-train]

# Jags model

model_code <- "
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = b1[g1[i]] + b2[g2[i]] # Note that b2 here is a parameter but b1 is not

    lh[i] <- dnorm(y[i], mu[i], sigma^-2)

    # Clustering model
    g2[i] ~ dcat(pi[1:B2])

    # Continuous environmental variable
    var2_1[i] ~ dnorm(mu_var2_1[g2[i]], sigma_var2_1^-2)
  }

  for(i in 1:N_test){
    y_test[i] ~ dnorm(mu_test[i], sigma^-2)
    mu_test[i] = b1[g1_test[i]] + b2[g2_test[i]]

    # Clustering model
    g2_test[i] ~ dcat(pi[1:B2])

    # Continuous variable
    var2_1_test[i] ~ dnorm(mu_var2_1[g2_test[i]], sigma_var2_1^-2)
   }

  # Priors

  # Prior on cluster membership
  pi ~ ddirch(alpha)

  for (g in 1:B2) {
    mu_var2_1_raw[g] ~ dnorm(0, 100^-2)
  }

  # Make sure these are in order to avoid label switching
  mu_var2_1 <- sort(mu_var2_1_raw[1:B2])

  # Prior on genotype effect
  for(i in 1:B1) {
    b1[i] ~ dnorm(0, sigma_b1^-2) # Prior on b1 effect
  }

  for(i in 1:B2) {
    b2[i] ~ dnorm(0, sigma_b2^-2) # Prior on b2 effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_var2_1 ~ dt(0, 10^-2, 1)T(0,)
  sigma_b1 ~ dt(0, 10^-2, 1)T(0,)
  sigma_b2 ~ dt(0, 10^-2, 1)T(0,)

  lp <- sum(lh) + sum(b1) + sum(b2)
  }
"


# Run model
usedB2 <- 3
model_data <- list(
  N = length(train),
  y = dat$df$y[train],
  y_test = dat$df$y[test],
  B2 = usedB2,
  B1 = B1,
  g1 = dat$df$g1[train],
  var2_1 = dat$df$var2_1[train],
  alpha = rep(1, usedB2),
  N_test = length(test),
  var2_1_test = dat$df$var2_1[test],
  g1_test = dat$df$g1[test]
)
str(model_data)

# Parameters to watch
model_parameters <- c("b1", "b2", "g2", "pi", "mu_var2_1", "mu", "sigma", "sigma_var2_1",
                      'g2_test','mu_test', 'lp')

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Results and output of the simulated example
plot(model_run)
print(model_run)

# RMSE train and test
caret::RMSE(dat$df$y[train], model_run$BUGSoutput$median$mu)
caret::RMSE(dat$df$y[test], model_run$BUGSoutput$median$mu_test)

# Plot the posterior cluster membership
qplot(dat$df$g2[train], model_run$BUGSoutput$median$g2) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_bw(base_size = 16) +
  labs(x = 'True cluster indicator', y = 'Estimated cluster indicator')

# Overall predictions
qplot(dat$df$y[test], model_run$BUGSoutput$median$mu_test) +
  geom_abline() +
  theme_bw(base_size = 16) +
  labs(x = 'y', y = expression(hat(y)))

# Prediction of b1 effects
qplot(dat$b1, model_run$BUGSoutput$median$b1) +
  geom_abline() +
  theme_bw(base_size = 16) +
  labs(x = expression(b[1]), y = expression(hat(b)[1])) +
# Prediction of b2 effects
qplot(dat$b2, model_run$BUGSoutput$median$b2) +
  geom_abline() +
  theme_bw(base_size = 16) +
  labs(x = expression(b[2]), y = expression(hat(b)[2]))




# Model 2: b1 + b2 + int --------------------------------------------------
# In this model we are clustering both variables


# Function to generate data
genData <- function(N, B1, B2, mu_var2_1, sigma, sigma_b1, sigma_b2, sigma_var2_1) {
  b1 <- rnorm(B1, 0, sigma_b1)
  b2 <- rnorm(B2, 0, sigma_b2)
  g1 <- rep(1:B1, each = N / B1)

  # cluster membership
  g2 <- sample(1:B2, size = N, replace = TRUE)

  # generate mean of the model
  mu_y <- b1[g1] + b2[g2]
  y <- rnorm(N, mu_y, sigma)
  var2_1 <- rnorm(N, mu_var2_1[g2], sigma_var2_1)
  df <- data.frame(y = y, g2 = g2, g1 = g1, var2_1 = var2_1)

  return(list(df = df, b1 = b1, b2 = b2))
}

# set up the values
N <- 400
B1 <- 10
B2 <- 3
mu_var2_1 <- c(-10, 0, 10)
sigma <- 1
sigma_b1 <- 1
sigma_b2 <- 1
sigma_var2_1 <- 1

# simulate a data set
set.seed(02)
dat <- genData(
  N = N, B1 = B1, B2 = B2, mu_var2_1 = mu_var2_1,
  sigma = sigma, sigma_b1 = sigma_b1,
  sigma_b2 = sigma_b2, sigma_var2_1 = sigma_var2_1
)

# plot data

hist(dat$df$var2_1, breaks = 30, freq = FALSE)
for (g in 1:B2) curve(dnorm(x, mean = mu_var2_1[g], sd = sigma_var2_1) / B2, col = g, add = TRUE)

# indexes for train and test
train <- sample(1:N, size = 200)
test <- (1:N)[-train]

# Jags model

model_code <- "
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = b1[g1[i]] + b2[g2[i]] # Note that b2 here is a parameter but b1 is not

    # Clustering model
    g2[i] ~ dcat(pi[1:B2])

    # Continuous environmental variable
    var2_1[i] ~ dnorm(mu_var2_1[g2[i]], sigma_var2_1^-2)
  }

  for(i in 1:N_test){
    y_test[i] ~ dnorm(mu_test[i], sigma^-2)
    mu_test[i] = b1[g1_test[i]] + b2[g2_test[i]]

    # Clustering model
    g2_test[i] ~ dcat(pi[1:B2])

    # Continuous variable
    var2_1_test[i] ~ dnorm(mu_var2_1[g2_test[i]], sigma_var2_1^-2)
   }

  # Priors

  # Prior on cluster membership
  pi ~ ddirch(alpha)

  for (g in 1:B2) {
    mu_var2_1_raw[g] ~ dnorm(0, 100^-2)
  }

  # Make sure these are in order to avoid label switching
  mu_var2_1 <- sort(mu_var2_1_raw[1:B2])

  # Prior on genotype effect
  for(i in 1:B1) {
    b1[i] ~ dnorm(0, sigma_b1^-2) # Prior on b1 effect
  }

  for(i in 1:B2) {
    b2[i] ~ dnorm(0, sigma_b2^-2) # Prior on b2 effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_var2_1 ~ dt(0, 10^-2, 1)T(0,)
  sigma_b1 ~ dt(0, 10^-2, 1)T(0,)
  sigma_b2 ~ dt(0, 10^-2, 1)T(0,)
  }
"


# Run model
usedB2 <- 3
model_data <- list(
  N = length(train),
  y = dat$df$y[train],
  y_test = dat$df$y[test],
  B2 = usedB2,
  B1 = B1,
  g1 = dat$df$g1[train],
  var2_1 = dat$df$var2_1[train],
  alpha = rep(1, usedB2),
  N_test = length(test),
  var2_1_test = dat$df$var2_1[test],
  g1_test = dat$df$g1[test]
)
str(model_data)

# Parameters to watch
model_parameters <- c("b1", "b2", "g2", "pi", "mu_var2_1", "mu", "sigma", "sigma_var2_1",
                      'g2_test','mu_test')

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Results and output of the simulated example
plot(model_run)
print(model_run)

# RMSE train and test
caret::RMSE(dat$df$y[train], model_run$BUGSoutput$median$mu)
caret::RMSE(dat$df$y[test], model_run$BUGSoutput$median$mu_test)

# Plot the posterior cluster membership
qplot(dat$df$g2[train], model_run$BUGSoutput$median$g2) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_bw(base_size = 16) +
  labs(x = 'True cluster indicator', y = 'Estimated cluster indicator')

# Overall predictions
qplot(dat$df$y[test], model_run$BUGSoutput$median$mu_test) +
  geom_abline() +
  theme_bw(base_size = 16) +
  labs(x = 'y', y = expression(hat(y)))

# Prediction of b1 effects
qplot(dat$b1, model_run$BUGSoutput$median$b1) +
  geom_abline() +
  theme_bw(base_size = 16) +
  labs(x = expression(b[1]), y = expression(hat(b)[1])) +
  # Prediction of b2 effects
  qplot(dat$b2, model_run$BUGSoutput$median$b2) +
  geom_abline() +
  theme_bw(base_size = 16) +
  labs(x = expression(b[2]), y = expression(hat(b)[2]))


# Model 1


samples_list <- model_run$BUGSoutput$sims.list
samples_list |> names()

post_class_p <- samples_list$g2

m = 3000
K = 3
J = 200

# initialize mcmc arrays
mcmc <- array(data = NA, dim = c(m = m, K = K, J = J))
mcmc |> str()

# assign posterior draws to the array
mcmc[, , 1] <- samples_list$pi
for (i in 1:(J - 1)) {
  mcmc[, , i + 1] <- runif(1, 0,1)
}

# set of selected relabeling algorithm
set <-
  c("PRA",
    "ECR",
    "ECR-ITERATIVE-1",
    "AIC",
    "ECR-ITERATIVE-2",
    "STEPHENS",
    "DATA-BASED")

set <-
  c("DATA-BASED")

# find the MAP draw as a pivot
mapindex = which.max(samples_list$lp)

ls_lcm <-
  label.switching(
    method = set,
    # zpivot = post_class_p[mapindex,],
    z = post_class_p,
    K = K,
    # prapivot = mcmc[mapindex, ,],
    # constraint = 1,
    # mcmc = mcmc,
    # p = post_class_p,
    data = dat$df$y[train]
  )

# ls_lcm <-
#   label.switching(
#     method = set,
#     zpivot = post_class_p[mapindex,],
#     z = post_class_p,
#     K = K,
#     prapivot = mcmc[mapindex, ,],
#     constraint = 1,
#     mcmc = mcmc,
#     p = post_class_p,
#     data = dat$df$y[train]
#   )

mcmc_permuted <- permute.mcmc(mcmc, ls_lcm$permutations$`DATA-BASED`)

mcmc_permuted$output[,,1][,1] |> plot()

