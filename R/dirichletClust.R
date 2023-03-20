# Dirichlet Process Clustering

library(R2jags)
library(ggplot2)
library(mclust)

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


G <- 5
N <- 500
muT <- c(-10,0,10,20,40) #runif(15, -10, 40)
sigma <- 1
se <- 1
st <- 1
I <- 5

set.seed(2022)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, se = se, st = st)
qplot(dat$df$t, geom = 'histogram', fill = I("steelblue"), colour = I("steelblue"), alpha = I(0.6)) + theme_bw() + labs(x = 'var env')
BIC <- mclustBIC(dat$df[,-2])
plot(BIC)

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
    mu_env[h] ~ dnorm(0, 10^-2)
    st[h] ~ dgamma(0.1, 0.1)
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

  for(i in 1:H) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on environment effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$df$y, t = dat$df$t, H = 9, a = 1)

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



# one numerical and one categorical variable

# generate data -----------------------------------------------------------

genData <- function(N, G, I, muT, sigma, sg, se, st){

  g <- rnorm(I, 0, sg)
  e <- rnorm(G, 0, se)
  gen <- rep(1:I, each = N/I)

  # cluster membership
  env <- sample(1:G, size = N, replace = TRUE)

  # generate mean of the model
  mu_y <- g[gen] + e[env]
  y <- rnorm(N, mu_y, sigma)
  t <- c(rnorm(N/2, muT[env], st), rnorm(N/2, muT[env], st))
  df <- data.frame(y = y, env = env, gen = gen, t = t)

  return(list(df = df, g = g, e = e))

}


G <- 15
N <- 3000
muT <- runif(15, -20, 50)
sigma <- 1
sg <- 1
se <- 1
st <- 1
I <- 5

set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st)
qplot(dat$df$t, geom = 'histogram', fill = I("steelblue"), colour = I("steelblue"), alpha = I(0.6)) + theme_bw() + labs(x = 'var env')


# Clustering function -----------------------------------------------------

model_code <- '
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] # Note that env here is a parameter but gen is not

    # Continuous environmental variable
    t[i] ~ dnorm(mu_env[env[i]], st[env[i]]^-2)
    # Clustering model
    env[i] ~ dcat(pi[])
   }

  # Priors
   for (h in 1:H) {
    mu_env[h] ~ dnorm(0, 10^-2)
    st[h] ~ dgamma(0.1, 0.1)
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

  # Prior on genotype effect
  for(i in 1:I) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:H) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on environment effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$df$y, I = I, gen = dat$df$gen,
                   t = dat$df$t, H = 20, a = 1)

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

# Create a plot of the posterior density line
post <- model_run$BUGSoutput
st <- post$mean$st
mu_env <- post$mean$mu_env
pi <- post$mean$pi

N <- length(dat$df$t)
df <- data.frame(t = dat$df$t)

df <- df %>%
  mutate(mixt = c(
    rnorm(n = N/2, mean = mu_env[1], sd = sqrt(st[1])),
    rnorm(n = N/2, mean = mu_env[2], sd = sqrt(st[2]))
  )) %>% as.data.frame()

df %>%
  ggplot(aes(t)) +
  geom_density(colour = "orange", size = 1) +
  geom_density(
    data = df,
    aes(mixt, y = ..density.. * pi[1]),
    colour = "royalblue", linetype = "dotted", size = 1.1
  ) +
  geom_density(
    data = df ,
    aes(mixt, y = ..density.. * pi[2]),
    colour = "plum", linetype = "dotted", size = 1.1
  ) +
  xlim(min(df$t) - 2, max(df$t) + 2) +
  labs(x = 'env comp', y = '') +
  # annotate("text",
  #          x = muT[1], y = 0.21, label = expression(mu[1]),
  #          size = 7
  # ) +
  # annotate("text",
  #          x = muT[2], y = 0.21, label = expression(mu[2]),
  #          size = 7
  # ) +
  theme_bw()




# # two numericals variables -----------------------------------------
# finalizar
rm(list = ls())

# one numerical and one categorical variable

# generate data -----------------------------------------------------------

genData <- function(N, G, I, muT, sigma, sg, se, st){

  g <- rnorm(I, 0, sg)
  e <- rnorm(G, 0, se)
  gen <- rep(1:I, each = N/I)

  # cluster membership
  env <- sample(1:G, size = N, replace = TRUE)

  # generate mean of the model
  mu_y <- g[gen] + e[env]
  y <- rnorm(N, mu_y, sigma)
  t <- c(rnorm(N/2, muT[env], st), rnorm(N/2, muT[env], st))
  df <- data.frame(y = y, env = env, gen = gen, t = t)

  return(list(df = df, g = g, e = e))

}


G <- 2
N <- 100
muT <- c(-5, 10)
sigma <- 1
sg <- 1
se <- 1
st <- 1
I <- 5

set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st)
qplot(dat$df$t, geom = 'histogram', fill = I("steelblue"), colour = I("steelblue"), alpha = I(0.6)) + theme_bw() + labs(x = 'var env')


# Clustering function -----------------------------------------------------

model_code <- '
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] # Note that env here is a parameter but gen is not

    # Continuous environmental variable
    t[i] ~ dnorm(mu_env[env[i]], st[env[i]]^-2)
    # Clustering model
    env[i] ~ dcat(pi[])
   }

  # Priors
   for (h in 1:H) {
    mu_env[h] ~ dnorm(0, 10^-2)
    st[h] ~ dgamma(0.1, 0.1)
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

  # Prior on genotype effect
  for(i in 1:I) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:H) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on environment effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$df$y, I = I, gen = dat$df$gen,
                   t = dat$df$t, H = 10, a = 1)

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
