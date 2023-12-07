# Clustering + Model
library(ggplot2)
library(R2jags)

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
  t <- rnorm(N, muT[env], st)
  df <- data.frame(y = y, env = env, gen = gen, t = t, blin = blin)

  return(list(df = df, g = g, e = e, Q = Q))

}


G <- 3
N <- 100
#muT <- c(-5, 10)
#muT <- c(-20,-5,5,10, 50)
set.seed(02)
muT <- runif(G,-50,100)
sigma <- 1
sg <- 10
se <- 10
st <- 1
I <- 10
lambda <- 12
set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st, lambda = lambda)

hist(dat$df$t, breaks = 30, freq = FALSE)
for (g in 1:G) curve(dnorm(x, mean = muT[g], sd = st)/G, col = g, add = TRUE)

# ggplot(dat$df, aes(x = t)) +
#   geom_histogram(aes(y = after_stat(density)),
#                  colour = 1, fill = "white") +
#   stat_function(fun = function(x) dnorm(x, mean = muT[1])/G, color = "steelblue", linewidth = 0.7) +
#   stat_function(fun = function(x) dnorm(x, mean = muT[2])/G, color = "firebrick", linewidth = 0.7) +
#   theme_bw() +
#   xlab('var')


## clustering

model_code_groups <- "
model
{
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu_g[Z[i]], sigma^-2)
    Z[i] ~ dcat(pi[i, 1:G])
    for (g in 1:G) {
      exp_theta[i, g] <- exp(theta[i, g])
      pi[i, g] <- exp(theta[i, g]) / sum(exp_theta[i, 1:G])
      theta[i, g] ~ dnorm(0, 6^-2)
    }
  }
  # Priors
  sigma ~ dt(0, 10^-2, 1)T(0,)
  for (g in 1:G) {
    mu_g_raw[g] ~ dnorm(0, 100^-2)
  }
  # Make sure these are in order to avoid label switching
  mu_g <- sort(mu_g_raw[1:G])
}
"

# Set up the data
usedG <- 10
model_data_groups <- list(N = N, y = dat$df$t, G = usedG)

# Choose the parameters to watch
model_parameters_groups <- c("mu_g", "sigma", "Z", "pi")

# Run the model
model_run_groups <- jags(
  data = model_data_groups,
  parameters.to.save = model_parameters_groups,
  model.file = textConnection(model_code_groups)
)

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
plot(model_run_groups)
print(model_run_groups)


model_run_groups$BUGSoutput$median$mu_g
Z <- model_run_groups$BUGSoutput$median$Z
as.integer(Z)

qplot(dat$df$env, model_run_groups$BUGSoutput$median$Z) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_light()


# Applying the model


model_code <- '
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] + blin[i] # Note that env here is a parameter but gen is not
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])
   }

  # Priors

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
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$df$y, I = I, gen = dat$df$gen, env = as.integer(Z), Q = dat$Q, G = usedG)

# Choose the parameters to watch
model_parameters <- c("g", "e", "mu", "sigma", "blin")

# Run the model
model_run_sep <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Results and output of the simulated example
plot(model_run)
print(model_run)

qplot(dat$df$y, model_run$BUGSoutput$median$mu) +
  geom_abline()+
  theme_light() +
  labs(x = expression(hat(y)), y = 'y')
round(caret::RMSE(dat$df$y, model_run_sep$BUGSoutput$median$mu), 4)



# generate the test data

G <- 10
N <- 3000
#muT <- c(-5, 10)
#muT <- c(-20,-5,5,10, 50)
set.seed(02)
muT <- runif(G,-50,100)
sigma <- 1
sg <- 10
se <- 10
st <- 1
I <- 10
lambda <- 12
set.seed(04)
dat_test <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st, lambda = lambda)
round(caret::RMSE(dat_test$df$y, model_run_sep$BUGSoutput$median$mu), 4)

RMSE_CplusBammit <- list()
for(i in c(1:100)){
  set.seed(i)
  dat_test <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st, lambda = lambda)
  RMSE_CplusBammit$N3000[i] <- round(caret::RMSE(dat_test$df$y, model_run_sep$BUGSoutput$median$mu), 4)
}

saveRDS(RMSE_CplusBammit, file = 'RMSE_CplusBammit')

