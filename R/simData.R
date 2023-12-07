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
N <- 500
muT <- c(-2,5, 10)
#muT <- c(-100,-50, -10, 0, 20, 50, 100, 200, 400, 500)
#muT <- c(3, 9)
sigma <- 1
sg <- 1
se <- 1
st <- 1
I <- 10

set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st)


str(dat)
table(dat$df$env)
hist(dat$df$t, breaks = 30, freq = FALSE)
for (g in 1:G) curve(dnorm(x, mean = muT[g], sd = st)/G, col = g, add = TRUE)

ggplot(dat$df, aes(x = t)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  stat_function(fun = function(x) dnorm(x, mean = muT[1])/G, color = "steelblue",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muT[2])/G, color = "firebrick",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muT[3])/G, color = "chartreuse4",size = 0.7) +
  theme_bw() +
  labs(x = 'temperature')

ggplot(dat$df, aes(x = y)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  theme_bw()


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

   for(i in 1:N_test){
   y_test[i] ~ dnorm(mu_test[i], sigma^-2)
    mu_test[i] = g[gen_test[i]] + e[env_test[i]] # Note that env here is a parameter but gen is not

    # Clustering model
    env_test[i] ~ dcat(pi[1:G])

    # Continuous environmental variable
    t_test[i] ~ dnorm(mu_env[env_test[i]], st^-2)
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
usedG <- 3
train <- sample(1:N, size = 200)
test <- (1:N)[-train]
model_data <- list(N = length(train), y = dat$df$y[train], G = usedG, I = I, gen = dat$df$gen[train],
                   t = dat$df$t[train], alpha = rep(1,usedG), N_test = length(test),
                   t_test = dat$df$t[test], gen_test = dat$df$gen[test], y_test = dat$df$y[test])
str(model_data)
# Choose the parameters to watch
model_parameters <- c("g", "e", "env", "pi", "mu_env", "mu", "sigma", "st", 'env_test',
                      'mu_test')

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

plot(model_run)
#traceplot(model_run, varname= 'env')

caret::RMSE(dat$df$y[train], model_run$BUGSoutput$median$mu)
caret::RMSE(dat$df$y[test], model_run$BUGSoutput$median$mu_test)




# Results and output of the simulated example
plot(model_run)
print(model_run)

# Plot the posterior cluster membership
qplot(dat$df$env[train], model_run$BUGSoutput$median$env) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_bw(base_size = 16) +
  labs(x = 'True cluster indicator', y = 'Estimated cluster indicator')

# Overall predictions
qplot(dat$df$y[train], model_run$BUGSoutput$median$mu) + geom_abline() +
  geom_linerange(aes(ymin =  apply(model_run$BUGSoutput$sims.list$mu, 2, function(x) quantile(x, 0.05)),
                     ymax = apply(model_run$BUGSoutput$sims.list$mu, 2, function(x) quantile(x, 0.95))), alpha = 0.5, size = 0.4) +
  theme_bw() +
  labs(x = 'True', y = 'Estimated')

qplot(dat$df$y[test], model_run$BUGSoutput$median$mu_test) +
  geom_abline() +
  geom_linerange(aes(ymin =  apply(model_run$BUGSoutput$sims.list$mu_test, 2, function(x) quantile(x, 0.05)),
                     ymax = apply(model_run$BUGSoutput$sims.list$mu_test, 2, function(x) quantile(x, 0.95))), alpha = 0.5, size = 0.4) +
  theme_bw() +
  labs(x = 'True', y = 'Estimated')


# Prediction of genotype effects
qplot(model_run$BUGSoutput$median$g, dat$g) + geom_abline()

# Prediction of environment effects
qplot(model_run$BUGSoutput$median$e, dat$e) + geom_abline()






