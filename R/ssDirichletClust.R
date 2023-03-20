# Spike and Slab Dirichlet Process Clustering

library(R2jags)
library(ggplot2)

rm(list = ls())


# just one variable

# generate data -----------------------------------------------------------

genData <- function(N, G, muT, sigma, se, st){

  e <- rnorm(G, 0, se)
  env <- sample(1:G, size = N, replace = TRUE)

  # generate mean of the model
  y <- rnorm(N, e[env], sigma)
  t <- c(rnorm(N/2, muT[env], st), rnorm(N/2, muT[env], st))
  df <- data.frame(y = y, env = env, t = t)

  return(list(df = df, e = e))

}


G <- 3
N <- 100
muT <- runif(3, -10, 40)
sigma <- 1
se <- 1
st <- 1
I <- 5

set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, se = se, st = st)
qplot(dat$df$t, geom = 'histogram', fill = I("steelblue"), colour = I("steelblue"), alpha = I(0.6)) + theme_bw() + labs(x = 'var env')


# Clustering function -----------------------------------------------------

model_code <- '
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(e[env[i]] , sigma^-2)

    # Continuous environmental variable
    t[i] ~ dnorm(mu_env[env[i]], st[env[i]]^-2)
    # Clustering model
    env[i] ~ dcat(pi[])
   }

  # Priors
   for (h in 1:H) {
    mu_env[h] ~ dnorm(mu_env_aux[h], tau_env)
    st[h] ~ dgamma(0.1, 0.1)
    mu_env_aux[h] = alpha + sum(beta[h]*V[h])
   }

  # Make sure these are in order to avoid label switching
  # mu_env <- sort(mu_env_raw[1:H])
  # st <- sort(st_raw[1:H])

  # Stick breaking prior
  for (h in 1:(H-1)) {
    V[h] ~ dbeta(1,a)
  }
  V[H] <- 1
  pi[1] <- V[1]
  for (h in 2:H) {
    pi[h] <- V[h] * (1-V[h-1]) * pi[h-1] / V[h-1]
  }

  for(j in 1:H){
    beta0[j] ~ dnorm(0,tau_env)
    beta1[j] ~ dbern(prob)       # Bernoulli distributed prior (so it can only be 1:included or 0:notincluded)
    beta[j] <- beta0[j]*beta1[j] # inclusion probability
  }
  prob ~ dunif(0,1)
  tau_env ~ dgamma(.1,.1)
  alpha ~ dnorm(0,0.1)

  for(i in 1:H) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on environment effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$df$y, t = dat$df$t, H = 10, a = 1)

# Choose the parameters to watch
model_parameters <- c("e", "env", "pi", "mu_env", "sigma", "st")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Results and output of the simulated example
plot(model_run)
