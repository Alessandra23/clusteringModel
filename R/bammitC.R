library(R2jags)
library(ggplot2)

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

genData <- function(N, G1, G2, G3, NG, NE, NV, muG, muE, muV, sg, se, sv, sigma, lambda){

  Q <- length(lambda)
  g <- rnorm(G1, 0, sg)
  e <- rnorm(G2, 0, se)
  v <- rnorm(G3, 0, sv)

  # g <- g - mean(g)
  # e <- e - mean(e)

  sgroupG <- MCMCpack::riwish(NG, diag(NG))
  sgroupE <- MCMCpack::riwish(NE, diag(NE))
  sgroupV <- MCMCpack::riwish(NV, diag(NV))

  # cluster gen
  gen <- rep(NA, N)
  g_group <- matrix(NA, ncol = NG, nrow = N)
  for (i in 1:N) {
    gen[i] <- sample(1:G1, size = 1, replace = TRUE)
    g_group[i, ] <- mvtnorm::rmvnorm(1, muG[gen[i], ], sgroupG)
  }

  # cluster env
  env <- rep(NA, N)
  e_group <- matrix(NA, ncol = NE, nrow = N)
  for (i in 1:N) {
    env[i] <- sample(1:G2, size = 1, replace = TRUE)
    e_group[i, ] <- mvtnorm::rmvnorm(1, muE[env[i], ], sgroupE)
  }

  var3 <- rep(NA, N)
  var3_group <- matrix(NA, ncol = NV, nrow = N)
  for (i in 1:N) {
    var3[i] <- sample(1:G3, size = 1, replace = TRUE)
    var3_group[i, ] <- mvtnorm::rmvnorm(1, muV[var3[i], ], sgroupV)
  }

  # generate lambda, gamma, delta and kappa
  gamma <- generateBlin(G1, Q)
  delta <- generateBlin(G2, Q)
  var3_int <- generateBlin(G3, Q)

  # generate bilinear term
  blin <- rep(0, N)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k] * gamma[gen, k] * delta[env, k] * var3_int[var3, k]
  }

  # generate mean of the model
  mu_y <- g[gen] + e[env] + v[var3] + blin
  y <- rnorm(N, mu_y, sigma)

  df <- data.frame(y = y, env = env, gen = gen, var3, blin = blin)

  return(list(df = df, g = g, e = e, v = v, g_group = g_group, e_group = e_group, var3_group = var3_group, Q = Q))

}


N <- 400
NG <- 2
NE <- 2
NV <- 2
sg <- 1
se <- 1
sv <- 1
#sgroup <-  diag(2)
sigma <- 1
lambda <- 20

trueC <- 3
G1 <- trueC
G2 <- trueC
G3 <- trueC

set.seed(02)
muG <- matrix(runif(NG*G1, -50,50), ncol = NG, nrow = G1)
muE <- matrix(runif(NE*G2, -50,50), ncol = NE, nrow = G2)
muV <- matrix(runif(NV*G3, -50,50), ncol = NV, nrow = G3)

dat <- genData(N = N, G1 = G1, G2 = G2, G3 = G3, NG = NG, NE = NE, NV = NV,
               muG = muG, muE = muE, muV = muV, sg = sg, se = se, sv = sv,
               sigma = sigma, lambda = lambda)

train <- sample(1:N, size = 200)
test <- (1:N)[-train]
usedC <-  3


# CBAMMIT -----------------------------------------------------------------

model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] + v[var3[i]] + blin[i] #gen and env are parameters
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q] * var3_int[var3[i],1:Q])

    # Clustering model
    gen[i] ~ dcat(piG[i, 1:G1])
    env[i] ~ dcat(piE[i, 1:G2])
    var3[i] ~ dcat(piV[i, 1:G3])

    # Continuous variables
    g_group[i, 1:NG] ~ dmnorm(mu_gen[gen[i], 1:NG], sgroup_inv_gen)
    e_group[i, 1:NE] ~ dmnorm(mu_env[env[i], 1: NE], sgroup_inv_env)
    var3_group[i, 1:NV] ~ dmnorm(mu_var3[var3[i], 1: NV], sgroup_inv_var3)

    for (j in 1:G1) {
      exp_theta_gen[i, j] <- exp(theta_gen[i, j])
      piG[i, j] <- exp(theta_gen[i, j]) / sum(exp_theta_gen[i, 1:G1])
      theta_gen[i, j] ~ dnorm(0, 6^-2)
    }


    for (j in 1:G2) {
      exp_theta_env[i, j] <- exp(theta_env[i, j])
      piE[i, j] <- exp(theta_env[i, j]) / sum(exp_theta_env[i, 1:G2])
      theta_env[i, j] ~ dnorm(0, 6^-2)
    }


    for (j in 1:G3) {
      exp_theta_var3[i, j] <- exp(theta_var3[i, j])
      piV[i, j] <- exp(theta_var3[i, j]) / sum(exp_theta_var3[i, 1:G3])
      theta_var3[i, j] ~ dnorm(0, 6^-2)
    }

   }


   for (i in 1:N_test) {
    # Model for phenotype
    y_test[i] ~ dnorm(mu[i], sigma^-2)
    mu_test[i] = g[gen_test[i]] + e[env_test[i]] + v[var3_test[i]] + blin_test[i] #gen and env are parameters
    blin_test[i] = sum(lambda[1:Q] * gamma[gen_test[i],1:Q] * delta[env_test[i],1:Q]* var3_int[var3_test[i],1:Q])

    # Clustering model
    gen_test[i] ~ dcat(piG[i, 1:G1])
    env_test[i] ~ dcat(piE[i, 1:G2])
    var3_test[i] ~ dcat(piV[i, 1:G3])

    # Continuous variables
    g_group_test[i, 1:NG] ~ dmnorm(mu_gen[gen_test[i], 1:NG], sgroup_inv_gen)
    e_group_test[i, 1:NE] ~ dmnorm(mu_env[env_test[i], 1: NE], sgroup_inv_env)
    var3_group_test[i, 1:NE] ~ dmnorm(mu_var3[var3_test[i], 1: NV], sgroup_inv_var3)

    for (j in 1:G1) {
      exp_theta_gen_test[i, j] <- exp(theta_gen_test[i, j])
      piG_test[i, j] <- exp(theta_gen_test[i, j]) / sum(exp_theta_gen_test[i, 1:G1])
      theta_gen_test[i, j] ~ dnorm(0, 6^-2)
    }


    for (j in 1:G2) {
      exp_theta_env_test[i, j] <- exp(theta_env_test[i, j])
      piE_test[i, j] <- exp(theta_env_test[i, j]) / sum(exp_theta_env_test[i, 1:G2])
      theta_env_test[i, j] ~ dnorm(0, 6^-2)
    }


    for (j in 1:G3) {
      exp_theta_var3_test[i, j] <- exp(theta_var3_test[i, j])
      piV_test[i, j] <- exp(theta_var3_test[i, j]) / sum(exp_theta_var3_test[i, 1:G3])
      theta_var3_test[i, j] ~ dnorm(0, 6^-2)
    }

   }

  # Priors on means gen
  for (j in 1:G1) {
    for (m in 1:NG) {
      mu_gen_raw[j, m] ~ dnorm(0, 100^-2)
    }
    for (m in 2:NG) {
      mu_gen[j, m] <- mu_gen_raw[j, m]
    }
  }

   # Sort first dimension to avoid label switching
   mu_gen[1:G1, 1] <- sort(mu_gen_raw[1:G1, 1])


   # Priors on means env
  for (j in 1:G2) {
    for (m in 1:NE) {
      mu_env_raw[j, m] ~ dnorm(0, 100^-2)
    }
    for (m in 2:NE) {
      mu_env[j, m] <- mu_env_raw[j, m]
    }
  }

   # Sort first dimension to avoid label switching
   mu_env[1:G2, 1] <- sort(mu_env_raw[1:G2, 1])



   # Priors on means var3
  for (j in 1:G3) {
    for (m in 1:NV) {
      mu_var3_raw[j, m] ~ dnorm(0, 100^-2)
    }
    for (m in 2:NV) {
      mu_var3[j, m] <- mu_var3_raw[j, m]
    }
  }

   # Sort first dimension to avoid label switching
   mu_var3[1:G3, 1] <- sort(mu_var3_raw[1:G3, 1])

   # Prior on genotype effect
  for(i in 1:G1) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G2) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  for(i in 1:G3) {
    v[i] ~ dnorm(0, sigma_v^-2) # Prior on var3 effect
  }

  # Priors on gamma
  for(q in 1:Q){
    for(i in 1:G1){
      thetaG[i,q] ~ dnorm(0,1)
    }
    mG[q] = sum(thetaG[1:G1,q])/G1
    for(i in 1:G1){
    thetaGNew[i,q] = thetaG[i,q] - mG[q]
    }
    sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:G1,q]^2 + 0.000001)))
    for(i in 1:G1){
      gamma[i,q] = thetaGNew[i,q]*sqrtThetaG[q]
    }
  }

   # Priors on delta
   for(q in 1:Q){
    for(j in 1:G2){
      thetaD[j,q] ~ dnorm(0,1)
    }
    mD[q] = sum(thetaD[1:G2,q])/G2
    for(j in 1:G2){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:G2,q]^2+ 0.000001)))
    for(j in 1:G2){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
   }

  # Priors on var3_int
   for(q in 1:Q){
    for(j in 1:G3){
      thetavar3[j,q] ~ dnorm(0,1)
    }
    mvar3[q] = sum(thetavar3[1:G3,q])/G3
    for(j in 1:G3){
    thetavar3New[j,q] = thetavar3[j,q] - mvar3[q]
    }
    sqrtThetavar3[q] = sqrt(1/(sum(thetavar3New[1:G3,q]^2+ 0.000001)))
    for(j in 1:G3){
      var3_int[j,q] = thetavar3New[j,q]*sqrtThetavar3[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, 100^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  sgroup_inv_gen ~ dwish(sG, NG)
  sgroup_gen <- inverse(sgroup_inv_gen)

  sgroup_inv_env ~ dwish(sE, NE)
  sgroup_env <- inverse(sgroup_inv_env)

  sgroup_inv_var3 ~ dwish(sV, NV)
  sgroup_var3 <- inverse(sgroup_inv_var3)

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
  sigma_v ~ dt(0, 10^-2, 1)T(0,)
}
"


# Run CBAMMIT

# Set up the data
model_data_CBAMMIT <- list(N = length(train), y = dat$df$y[train], G1 = usedC, G2 = usedC, G3 = usedC,
                           NG = NG, NE = NE, NV = NV, sG = diag(NG), sE = diag(NE), sV = diag(NV), Q = dat$Q,
                           g_group = dat$g_group[train,], e_group = dat$e_group[train,], var3_group = dat$var3_group[train,],
                           N_test = length(test), g_group_test = dat$g_group[test,],
                           e_group_test = dat$e_group[test,], var3_group_test = dat$var3_group[test,]) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
str(model_data_CBAMMIT)

# Choose the parameters to watch
model_parameters <- c("g", "e", "gen", "env", "piG", "piE", "mu_env", "mu_gen", "mu_var3",
                      "mu", 'mu_test', 'blin_test', "sigma", "sgroup_gen", 'sgroup_env', 'blin', 'gamma', 'delta', 'lambda')

# Run the model
model_run <- R2jags::jags(
  data = model_data_CBAMMIT,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

plot(model_run)

caret::RMSE(dat$df$y[train], model_run$BUGSoutput$median$mu)
caret::RMSE(dat$df$y[test], model_run$BUGSoutput$median$mu_test)
