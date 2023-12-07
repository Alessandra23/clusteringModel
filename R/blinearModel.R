library(R2jags)
library(ggplot2)
library(patchwork)

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


G <- 2
N <- 100
muT <- c(3, 9)
#muT <- c(30, 50)
set.seed(02)
muT <- runif(20,-10,30)
sigma <- 1
sg <- 1
se <- 1
st <- 1
I <- 10
lambda <- 100
set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, I = I, sg = sg, se = se, st = st, lambda = lambda)

hist(dat$df$t, breaks = 30, freq = FALSE)
for (g in 1:G) curve(dnorm(x, mean = muT[g], sd = st)/G, col = g, add = TRUE)

ggplot(dat$df, aes(x = t)) +
  geom_histogram(aes(y = after_stat(density)),
                 colour = 1, fill = "white") +
  stat_function(fun = function(x) dnorm(x, mean = muT[1])/G, color = "steelblue", linewidth = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muT[2])/G, color = "firebrick", linewidth = 0.7) +
  theme_bw() +
  xlab('Temperature')


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
Gused = 2
# Set up the data
model_data <- list(N = N, y = dat$df$y, G = Gused, I = I, gen = dat$df$gen, Q = dat$Q,
                   t = dat$df$t, alpha = rep(1,Gused))

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

# Plot the posterior cluster membership
qplot(dat$df$env, model_run$BUGSoutput$median$env) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_light() +
  labs(y = "Estimated environment cluster", x = 'True environment cluster')


# Overall predictions
p1 <- qplot(dat$df$y, model_run$BUGSoutput$median$mu) +
  geom_abline()+
  theme_light() +
  labs(x = expression(hat(y)), y = 'y')
p1
caret::RMSE(dat$df$y, model_run$BUGSoutput$median$mu)

# Prediction of genotype effects
g_plot <- qplot(dat$g, model_run$BUGSoutput$median$g) +
  geom_abline() + theme_bw() +
  labs(y = expression(hat(g)), x = 'g')

# Prediction of environment effects
e_plot <- qplot(dat$e, model_run$BUGSoutput$median$e) +
  geom_abline() + theme_bw() +
  labs(y = expression(hat(e)), x = 'e')

caret::RMSE(model_run$BUGSoutput$median$e, dat$e)

# Prediction of interactio effect
p2 <- qplot(dat$df$blin, model_run$BUGSoutput$median$blin) +
  geom_abline() + theme_light() +
  labs(y = "Estimated interaction", x = 'True interaction')
p2

p1+p2

g_plot + e_plot + p2
