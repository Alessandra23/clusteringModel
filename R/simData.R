library(R2jags)
library(ggplot2)

rm(list = ls())

genData <- function(N, G, I, muT, sigma, sg, se, st){

  g <- rnorm(I, 0, sg)
  e <- rnorm(G, 0, se)
  gen <- rep(1:I, each = N/I)

  # cluster membership
  env <- sample(1:G, size = N, replace = TRUE)

  # generate mean of the model
  mu_y <- g[gen] + e[env]
  y <- rnorm(N, mu_y, sigma)
  t <- rnorm(N, muT[env], st)
  df <- data.frame(y = y, env = env, gen = gen, t = t)

  return(list(df = df, g = g, e = e))

}


G <- 3
N <- 100
muT <- c(-5, 10, 30)
sigma <- 1
sg <- 10
se <- 10
st <- 1
I <- 5

set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st)

hist(dat$df$t, breaks = 30, freq = FALSE)
for (g in 1:G) curve(dnorm(x, mean = muT[g], sd = st)/G, col = g, add = TRUE)

# ggplot(dat, aes(x = t)) +
#   geom_histogram(aes(y = ..density..),
#                  colour = 1, fill = "white") +
#   stat_function(fun = function(x) dnorm(x, mean = muT[1])/G, color = "steelblue",size = 0.7) +
#   stat_function(fun = function(x) dnorm(x, mean = muT[2])/G, color = "firebrick",size = 0.7) +
#   stat_function(fun = function(x) dnorm(x, mean = muT[3])/G, color = "chartreuse4",size = 0.7) +
#   theme_bw()
#
# ggplot(dat, aes(x = y)) +
#   geom_histogram(aes(y = ..density..),
#                  colour = 1, fill = "white") +
#   theme_bw()


# Clustering function -----------------------------------------------------


model_code <- '
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] # Note that env here is a parameter but gen is not

    # Clustering model
    env[i] ~ dcat(pi[1:G])

    # Continuous environmental variable
    t[i] ~ dnorm(mu_env[env[i]], st^-2)
   }

  # Priors

  # Prior on cluster membership
  pi ~ ddirch(alpha)

  for (g in 1:G) {
    mu_env_raw[g] ~ dnorm(0, 100^-2)
  }
  # Make sure these are in order to avoid label switching
  mu_env <- sort(mu_env_raw[1:G])

  # Prior on genotype effect
  for(i in 1:I) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  st ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$df$y, G = G, I = I, gen = dat$df$gen,
                   t = dat$df$t, alpha = rep(1,G))

# Choose the parameters to watch
model_parameters <- c("g", "e", "env", "pi", "mu_env", "mu", "sigma", "st")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Results and output of the simulated example
plot(model_run)
print(model_run)

# Plot the posterior cluster membership
qplot(model_run$BUGSoutput$median$env, dat$df$env) +
  geom_jitter(width = 0.1, height = 0.1)

# Overall predictions
qplot(model_run$BUGSoutput$median$mu, dat$df$y) + geom_abline()

# Prediction of genotype effects
qplot(model_run$BUGSoutput$median$g, dat$g) + geom_abline()

# Prediction of environment effects
qplot(model_run$BUGSoutput$median$e, dat$e) + geom_abline()
