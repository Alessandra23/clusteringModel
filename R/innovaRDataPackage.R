#devtools::install_github("danilosarti/InnoVaR")
library(InnoVaR)
library(ggplot2)
library(tidyr)
library(R2jags)

rm(list = ls())

simInnoVar <- readRDS("~/Documents/GitHub/bammit/Real data/simInnoVar.RDS")

data <- simInnoVar[,c('environment', 'genotype', 'P', 'K', 'N', 'Clay', 'target')]
colnames(data)[7] <- "y"


set.seed(2022)
sample_size <- 80
sample_meanvector <- c(5.3, 30)
sample_covariance_matrix <- matrix(c(10, 6, 4, 10),
                                   ncol = 2)
# create bivariate normal distribution
sample_distribution1 <- MASS::mvrnorm(n = sample_size,
                               mu = sample_meanvector,
                               Sigma = sample_covariance_matrix)

sample_distribution2 <- MASS::mvrnorm(n = sample_size,
                                      mu = sample_meanvector,
                                      Sigma = sample_covariance_matrix)
cor(sample_distribution1)

data$y <- c(sample_distribution1[,1], sample_distribution2[,1])
data$Clay <- c(sample_distribution1[,2], sample_distribution2[,2])

cor(data[,-c(1,2)])

# histograms of the numeric variables
ggplot(gather(data[,-c(1,2)]), aes(value)) +
  geom_histogram(bins = 10) +
  facet_wrap(~key, scales = 'free_x') +
  theme_bw()

levels(data$genotype)


# One numerical

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
N <- length(data$y)
G <- 2
y <- data$y
I <- length(unique(data$genotype))
gen <- data$genotype
t <- data$Clay


model_data <- list(N = N, y = y, G = G, I = I, gen = gen,
                   t = t, alpha = rep(1,G))

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

# bilinear model


model_code <- '
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] + blin[i] # Note that env here is a parameter but gen is not
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

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

  sigma ~ dt(0, 10^-2, 1)T(0,)
  st ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data

N <- length(data$y)
G <- 2
y <- data$y
I <- length(unique(data$genotype))
gen <- data$genotype
t <- data$Clay
Q <-  1

model_data <- list(N = N, y = y, G = G, I = I, gen = gen, Q = Q,
                   t = t, alpha = rep(1,G))

# Choose the parameters to watch
model_parameters <- c("g", "e", "env", "pi", "mu_env", "mu", "sigma", "st", "blin")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Results and output of the simulated example
plot(model_run)
print(model_run)

model_run$BUGSoutput$mean$env


# Plot the posterior cluster membership
qplot(model_run$BUGSoutput$median$env, data$environment) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_light() +
  labs(x = "Estimated environment cluster", y = 'True environment cluster')




# Overall predictions
qplot(model_run$BUGSoutput$median$mu, y) +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_abline()+
  theme_light() +
  labs(x = expression(hat(y)), y = 'y')


# Prediction of genotype effects
qplot(model_run$BUGSoutput$median$g, gen) + geom_abline()

# Prediction of environment effects
qplot(model_run$BUGSoutput$median$e, data$environment) + geom_abline()


# multivariate model


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

N <- length(data$y)
G <- 2
NE <- 2
y <- data$y
I <- length(unique(data$genotype))
gen <- data$genotype
t <- data[, c('Clay', 'P')]#, 'K', 'N')]
Q <-  1

model_data <- list(N = N, y = y, G = G, I = I, NE = NE, gen = gen,
                   t = t, stE = diag(NE), Q = Q)#, alpha = rep(1,G))

# Choose the parameters to watch
model_parameters <- c("g", "e", "env", "pi", "muT", "mu", "sigma", "st", 'blin')
#model_parameters <- c("g")

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
qplot(model_run$BUGSoutput$median$env, env) +
  geom_jitter(width = 0.1, height = 0.1) + theme_light() +
  labs(x = 'environment cluster', y = 'environment cluster estimated')

# Overall predictions
qplot(model_run$BUGSoutput$median$mu, y) +
  geom_abline() + theme_light() +
  labs(x = expression(mu), y = 'y')





