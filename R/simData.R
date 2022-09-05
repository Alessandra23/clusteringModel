library(R2jags)
library(ggplot2)

rm(list = ls())

genData <- function(N, G, I, muT, sigma, sg,se, st){

  gen <- rnorm(I, 0, sg)
  env <- rnorm(G, 0, se)
  g <- rep(1:I, each = N/I)

  # clustering Temp
  theta <- matrix(rnorm(N * G, 0, 3), ncol = G, nrow = N)
  pi <- exp(theta) / apply(exp(theta), 1, sum)
  Z <- rep(NA, N)
  for (i in 1:N) Z[i] <- sample(1:G, size = 1, prob = pi[i, ])

  # generate mean of the model
  mu_y <- gen[g] + env[Z]
  y <- rnorm(N, mu_y, sigma)
  t <- rnorm(N,muT[Z], st)
  df <- data.frame(y = y, Z = Z, g = g, t = t)

  return(df)
  #blin = blin,
  #Q = Q))

}


G <- 3
N <- 100
muT <- c(-5, 10, 30)
sigma <- 10
sg <-  100
se = 0.001
st = 1
I <-  5
#lambda <- 12
set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st)
hist(dat$t, breaks = 30, freq = FALSE)
for (g in 1:G) curve(dnorm(x, mean = muT[g])/G, col = g, add = TRUE)

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
  model
  {
  # Likelihood
   for (i in 1:N) {
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = gen[genotype[i]] + env[Z[i]]
    Z[i] ~ dcat(pi[1:G])
    t[i] ~ dnorm(muT[Z[i]], st^-2)
   }

  # Priors
  pi ~ ddirch(alpha)

  for (g in 1:G) {
    muT_raw[g] ~ dnorm(0, 100^-2)
  }

  # Make sure these are in order to avoid label switching
  muT <- sort(muT_raw[1:G])


  # Prior on genotype effect
  for(i in 1:I) {
  gen[i] ~ dnorm(0, 100^-2) # Prior on genotype effect
  }

  #
  for(i in 1:G) {
  env[i] ~ dnorm(0, 100^-2) # Prior on genotype effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  st ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$y, G = G, I = I, genotype = dat$g, t = dat$t, alpha = rep(1,G))

# Choose the parameters to watch
model_parameters <- c("muT", "sigma", "Z", "pi", "mu", "env", "gen", "st")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
plot(model_run)
print(model_run)

qplot(model_run$BUGSoutput$median$Z, dat$Z)
qplot(model_run$BUGSoutput$median$mu, dat$y)
boxplot(dat$y~dat$Z)
boxplot(dat$y~dat$g)











