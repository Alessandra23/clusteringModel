# C + AMBARTI -----------------------------------------------------------------
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

genData <- function(N, G1, G2, NG, NE, muG, muE, sg, se, sigma, lambda){

  Q <- length(lambda)
  g <- rnorm(G1, 0, sg)
  e <- rnorm(G2, 0, se)

  sgroupG <- MCMCpack::riwish(NG, diag(NG))
  sgroupE <- MCMCpack::riwish(NE, diag(NE))

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

  # generate lambda, gamma, delta and kappa
  gamma <- generateBlin(G1, Q)
  delta <- generateBlin(G2, Q)

  # generate bilinear term
  blin <- rep(0, N)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k] * gamma[gen, k] * delta[env, k]
  }

  # generate mean of the model
  mu_y <- g[gen] + e[env] + blin
  y <- rnorm(N, mu_y, sigma)

  df <- data.frame(y = y, env = env, gen = gen, blin = blin)

  return(list(df = df, g = g, e = e, g_group = g_group, e_group = e_group, Q = Q))

}

simData <- function(N = 120, trueC, NE = 3, NG = 3, seedTrain = 02, seedTest = 04){

  G1 <- trueC
  G2 <- trueC
  N <- N
  NE <- NE
  NG <- NG
  sg <- 1
  se <- 1
  #sgroup <-  diag(2)
  sigma <- 1
  lambda <- 100

  set.seed(seedTrain)
  muG <- matrix(runif(NG*G1, -50,50), ncol = NG, nrow = G1)
  muE <- matrix(runif(NE*G2, -50,50), ncol = NE, nrow = G2)
  datTrain <- genData(N = N, G1 = G1, G2 = G2, NG = NG, NE = NE,
                      muG = muG, muE = muE, sg = sg, se = se,
                      sigma = sigma, lambda = lambda)

  set.seed(seedTest)
  muG <- matrix(runif(NG*G1, -50,50), ncol = NG, nrow = G1)
  muE <- matrix(runif(NE*G2, -50,50), ncol = NE, nrow = G2)
  datTest <- genData(N = N, G1 = G1, G2 = G2, NG = NG, NE = NE,
                     muG = muG, muE = muE, sg = sg, se = se,
                     sigma = sigma, lambda = lambda)

  return(list(datTrain = datTrain,
              datTest = datTest))

}

N <- 50
NG <- 2
NE <- 2
trueC <- 2
data <- simData(N = N, trueC = trueC, NE = NE, NG = NG, seedTrain = 02, seedTest = 04)
dat <- data$datTrain
usedC <-  6



# C + BAMMIT

model_code_clustG <- "
model
{
  # Likelihood
   for (i in 1:N) {
    gen[i] ~ dcat(piG[i, 1:G1])

    # Continuous variables
    g_group[i, 1:NG] ~ dmnorm(mu_gen[gen[i], 1:NG], sgroup_inv_gen)

    for (j in 1:G1) {
      exp_theta_gen[i, j] <- exp(theta_gen[i, j])
      piG[i, j] <- exp(theta_gen[i, j]) / sum(exp_theta_gen[i, 1:G1])
      theta_gen[i, j] ~ dnorm(0, 6^-2)
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

  sgroup_inv_gen ~ dwish(sG, NG)
  sgroup_gen <- inverse(sgroup_inv_gen)

}
"

model_code_clustE <- "
model
{
  # Likelihood
  for (i in 1:N) {
    env[i] ~ dcat(piE[i, 1:G2])

    # Continuous variables
    e_group[i, 1:NE] ~ dmnorm(mu_env[env[i], 1: NE], sgroup_inv_env)

    for (j in 1:G2) {
      exp_theta_env[i, j] <- exp(theta_env[i, j])
      piE[i, j] <- exp(theta_env[i, j]) / sum(exp_theta_env[i, 1:G2])
      theta_env[i, j] ~ dnorm(0, 6^-2)
    }

   }

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

  sgroup_inv_env ~ dwish(sE, NE)
  sgroup_env <- inverse(sgroup_inv_env)
}
"

model_data_groupsG <- list(N = N, G1 = usedC, NG = NG,
                           g_group = dat$g_group, sG = diag(NG))
model_data_groupsE <- list(N = N, G2 = usedC, NE = NE,
                           e_group = dat$e_group, sE = diag(NE))

# Choose the parameters to watch
model_parameters_groupsG <- c("gen", "piG","mu_gen", "sgroup_gen")
model_parameters_groupsE <- c("env", "piE", "mu_env", 'sgroup_env')

# Run the model
model_run_groupsG <- jags(
  data = model_data_groupsG,
  parameters.to.save = model_parameters_groupsG,
  model.file = textConnection(model_code_clustG)
)

model_run_groupsE <- jags(
  data = model_data_groupsE,
  parameters.to.save = model_parameters_groupsE,
  model.file = textConnection(model_code_clustE)
)


groupsG <- as.integer(model_run_groupsG$BUGSoutput$median$gen)
groupsE <- as.integer(model_run_groupsE$BUGSoutput$median$env)


# Run BAMMIT

model_code_BAMMIT <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = muall + g[gen[i]] + e[env[i]] + blin[i] #gen and env are parameters
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

   }

   muall ~ dnorm(mmu, smu^-2) # grand mean

   # Prior on genotype effect
  for(i in 1:G1) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G2) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
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

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, 100^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
}
"


# Set up the data
model_data_BAMMIT <- list(N = N, y = dat$df$y,
                          G1 = usedC, G2 = usedC,
                          gen = groupsG, env = groupsE,
                          Q = dat$Q, mmu = 100,
                          smu = 10) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
#str(model_data_BAMMIT)
# Choose the parameters to watch
model_parameters_BAMMIT <- c("g", "e", "mu", "sigma",  'blin', 'muall')

# Run the model
model_run_BAMMIT <- R2jags::jags(
  data = model_data_BAMMIT,
  parameters.to.save = model_parameters_BAMMIT,
  model.file = textConnection(model_code_BAMMIT)
)

round(caret::RMSE(model_run_BAMMIT$BUGSoutput$median$mu, dat$df$y), 4)
round(caret::RMSE(model_run_BAMMIT$BUGSoutput$median$mu, data$datTest$df$y), 4)


predFunc <- function(model, groupsG, groupsE){

  muhat <- model$mean$muall
  b1 <- model$mean$g
  b2 <- model$mean$e
  inthat <- model$mean$blin

  N <- length(groupsG)
  yhat <- rep(muhat, N) + b1[groupsG] + b2[groupsE] + inthat
  return(yhat)
}

predy <- predFunc(model = model_run_BAMMIT$BUGSoutput, groupsG = groupsG, groupsE = groupsE)
round(caret::RMSE(predy, data$datTest$df$y), 4)

# run ambarti
library(AMBARTI)

df = dat$df
df$g = groupsG
df$e = groupsE
df$y = df$y
df = df[,-which(colnames(df) %in% c('gen','env','blin'))]
df = as.data.frame(df)
y = df$y
#y <- rnorm(120,10,2)
x <-  df[,-which(colnames(df) == 'y')]
fit.ambarti.ammi = ambarti(x, y, ntrees = 50, nburn = 500, npost = 1000, nsteps = 1)
#saveRDS(fit.ambarti.ammi, "~/Documents/GitHub/bammit/Running models/fit.ambarti.ammi.RData")

qq = var_used_trees(fit.ambarti.ammi)
df2 = data$datTest$df
df2$g = gsub('g','', df2$gen)
df2$e = gsub('e','', df2$env)
df2$y = df2$y
df2 = df2[,-which(colnames(df2)%in% c('gen','env','blin'))]
df2 = as.data.frame(df2)
y_test = df2$y
x_test = df2[,-which(colnames(df2) == 'y')]
yhat_ambarti2_ammi = predict_ambarti_alessa(object=fit.ambarti.ammi, newdata = x_test, type = 'mean')
#saveRDS(yhat_ambarti2_ammi, "~/Documents/GitHub/bammit/Running models/yhat_ambarti2.RData")

caret::RMSE(y_test, yhat_ambarti2_ammi[,1])
caret::R2(y_test, yhat_ambarti2_ammi[,1])



# CBAMMIT -----------------------------------------------------------------

model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = g[gen[i]] + e[env[i]] + blin[i] #gen and env are parameters
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

    # Clustering model
    gen[i] ~ dcat(piG[i, 1:G1])
    env[i] ~ dcat(piE[i, 1:G2])

    # Continuous variables
    g_group[i, 1:NG] ~ dmnorm(mu_gen[gen[i], 1:NG], sgroup_inv_gen)
    e_group[i, 1:NE] ~ dmnorm(mu_env[env[i], 1: NE], sgroup_inv_env)

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

   # Prior on genotype effect
  for(i in 1:G1) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G2) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
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

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, 100^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  sgroup_inv_gen ~ dwish(sG, NG)
  sgroup_gen <- inverse(sgroup_inv_gen)

  sgroup_inv_env ~ dwish(sE, NE)
  sgroup_env <- inverse(sgroup_inv_env)

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
}
"


# Set up the data
model_data <- list(N = N, y = dat$df$y, G1 = usedC, G2 = usedC, NG = NG, NE = NE, gen = dat$df$gen, env = dat$df$env,
                   g_group = dat$g_group, e_group = dat$e_group, sG = diag(NG), sE = diag(NE), Q = dat$Q) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
str(model_data)
# Choose the parameters to watch
model_parameters <- c("g", "e", "gen", "env", "piG", "piE", "mu_env", "mu_gen",
                      "mu", "sigma", "sgroup_gen", 'sgroup_env', 'blin')

# Run the model
model_run <- R2jags::jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

model_run$BUGSoutput$median$gen

round(caret::RMSE(model_run$BUGSoutput$median$mu, dat$df$y), 4)
round(caret::RMSE(model_run$BUGSoutput$median$mu, data$datTest$df$y), 4)


