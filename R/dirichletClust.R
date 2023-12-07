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
muT <- c(-50,-10,0,50,100) #runif(15, -10, 40)
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


G <- 10
N <- 200
muT <- runif(10, -200, 500)
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

caret::RMSE(dat$df$y, model_run$BUGSoutput$median$mu)

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





# With interaction --------------------------------------------------------

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

# one numerical and one categorical variable

# generate data -----------------------------------------------------------

genData <- function(N, G, I, muT, sigma, sg, se, st, lambda){

  Q <- length(lambda)
  g <- rnorm(I, 0, sg)
  e <- rnorm(G, 0, se)
  gen <- rep(1:I, each = N/I)

  # cluster membership
  env <- sample(1:G, size = N, replace = TRUE)

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
  t <- c(rnorm(N/2, muT[env], st), rnorm(N/2, muT[env], st))
  df <- data.frame(y = y, env = env, gen = gen, t = t, blin = blin)

  return(list(df = df, g = g, e = e, Q = Q))

}


G <- 20
N <- 200
muT <- runif(G, -200, 500)
sigma <- 1
sg <- 1
se <- 1
st <- 1
I <- 5
lambda <- 10

set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st, lambda = lambda)
qplot(dat$df$t, geom = 'histogram', fill = I("steelblue"), colour = I("steelblue"), alpha = I(0.6)) + theme_bw() + labs(x = 'var env')
plot(dat$df$t)


train <- sample(1:N, size = 100)
test <- (1:N)[-train]


# Clustering function -----------------------------------------------------

model_code <- '
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] + blin[i] # Note that env here is a parameter but gen is not
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

    # Continuous environmental variable
    t[i] ~ dnorm(mu_env[env[i]], st[env[i]]^-2)
    # Clustering model
    env[i] ~ dcat(pi[])
  }


   for (i in 1:N_test) {
    # Model for phenotype
    y_test[i] ~ dnorm(mu_test[i], sigma^-2)
    mu_test[i] = g[gen_test[i]] + e[env_test[i]] + blin_test[i] # Note that env here is a parameter but gen is not
    blin_test[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

    # Continuous environmental variable
    t_test[i] ~ dnorm(mu_env[env_test[i]], st[env_test[i]]^-2)
    # Clustering model
    env_test[i] ~ dcat(pi[])
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
    for(j in 1:H){
      thetaD[j,q] ~ dnorm(0,1)
    }
    mD[q] = sum(thetaD[1:H,q])/H
    for(j in 1:H){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:H,q]^2+ 0.000001)))
    for(j in 1:H){
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
'

# Set up the data
model_data <- list(N = length(train), y = dat$df$y[train], I = I, gen = dat$df$gen[train], Q = dat$Q,
                   t = dat$df$t[train], H = 50, a = 1, N_test = length(test), t_test = dat$df$t[test],
                   gen_test = dat$df$gen[train])

# Choose the parameters to watch
model_parameters <- c("g", "e", "env", "pi", "mu_env", "mu", "sigma", "st", 'blin', 'mu_test')

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Results and output of the simulated example
plot(model_run)
#print(model_run)

caret::RMSE(dat$df$y[train], model_run$BUGSoutput$median$mu)
caret::RMSE(dat$df$y[test], model_run$BUGSoutput$median$mu_test)
