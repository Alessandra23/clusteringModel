# C + AMBARTI -----------------------------------------------------------------
library(R2jags)

# library(devtools)
# install_github("ebprado/AMBARTI/R package", ref='main')

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
  #browser()
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
trueC <- 3
data <- simData(N = N, trueC = trueC, NE = NE, NG = NG, seedTrain = 02, seedTest = 04)
dat <- data$datTrain
usedC <-  4

x_train <- dat$df
x_test <- data$datTest$df

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


    for (i in 1:N_test) {
    gen_test[i] ~ dcat(piG_test[i, 1:G1])

    # Continuous variables
    g_group_test[i, 1:NG] ~ dmnorm(mu_gen[gen_test[i], 1:NG], sgroup_inv_gen)

    for (j in 1:G1) {
      exp_theta_gen_test[i, j] <- exp(theta_gen_test[i, j])
      piG_test[i, j] <- exp(theta_gen_test[i, j]) / sum(exp_theta_gen_test[i, 1:G1])
      theta_gen_test[i, j] ~ dnorm(0, 6^-2)
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


   for (i in 1:N_test) {
    env_test[i] ~ dcat(piE_test[i, 1:G2])

    # Continuous variables
    e_group_test[i, 1:NE] ~ dmnorm(mu_env[env_test[i], 1: NE], sgroup_inv_env)

    for (j in 1:G2) {
      exp_theta_env_test[i, j] <- exp(theta_env_test[i, j])
      piE_test[i, j] <- exp(theta_env_test[i, j]) / sum(exp_theta_env_test[i, 1:G2])
      theta_env_test[i, j] ~ dnorm(0, 6^-2)
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
    mu[i] =  g[gen[i]] + e[env[i]] + blin[i] #gen and env are parameters
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

   }

   for (i in 1:N_test) {
    # Model for phenotype
    y_test[i] ~ dnorm(mu_test[i], sigma^-2)
    mu_test[i] =  g[gen_test[i]] + e[env_test[i]] + blin_test[i] #gen and env are parameters
    blin_test[i] = sum(lambda[1:Q] * gamma[gen_test[i],1:Q] * delta[env_test[i],1:Q])

   }

   #muall ~ dnorm(mmu, smu^-2) # grand mean

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

# groupsG <- dat$df$gen
# groupsE <- dat$df$env

# Set up the data
model_data_BAMMIT <- list(N = N, y = dat$df$y,
                          G1 = usedC, G2 = usedC,
                          gen = groupsG, env = groupsE,
                          Q = dat$Q)#, mmu = 100,
                          #smu = 10) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
#
# model_data_BAMMIT <- list(N = N, y = dat$df$y,
#                           G1 = usedC, G2 = usedC,
#                           gen = dat$df$gen, env = dat$df$env,
#                           Q = dat$Q) #, mmu = 100,
#                           #smu = 10)
str(model_data_BAMMIT)

# Choose the parameters to watch
model_parameters_BAMMIT <- c("g", "e", "sigma", 'mu', 'blin', 'lambda', 'gamma', 'delta')

# Run the model
model_run_BAMMIT <- R2jags::jags(
  data = model_data_BAMMIT,
  parameters.to.save = model_parameters_BAMMIT,
  model.file = textConnection(model_code_BAMMIT)
)


#
x_train
y_train <- data.frame(y = model_run_BAMMIT$BUGSoutput$median$mu,
                      g = groupsG,
                      e = groupsE)


#plot(model_run_BAMMIT)

round(caret::RMSE(model_run_BAMMIT$BUGSoutput$median$mu, dat$df$y), 4)
round(caret::RMSE(model_run_BAMMIT$BUGSoutput$median$mu, data$datTest$df$y), 4)


predFunc <- function(model = model_run_BAMMIT$BUGSoutput,
                     indexTrain = list(groupsG = groupsG,
                                       groupsE = groupsE),
                     indexTest = list(gtest = data$datTest$df$gen,
                                      etest = data$datTest$df$env)){

  #muhat <- model$mean$muall
  b1 <- model$mean$g
  b2 <- model$mean$e

  groupsG <-  indexTrain$groupsG
  groupsE <- indexTrain$groupsE
  aux <- cbind(gen = groupsG, env = groupsE, int = model_run_BAMMIT$BUGSoutput$mean$blin)
  inthat <- aux |> as.data.frame() |>
    dplyr::distinct(gen, env,.keep_all = TRUE)


  gtest <-  indexTest$gtest
  etest <- indexTest$etest

  aux2 <- cbind(gen = groupsG, env = groupsE) |> as.data.frame()
  int <- plyr::join(aux2, inthat)
  int[is.na(int)] <- 0


  N <- length(groupsG)
  yhat <- b1[groupsG] + b2[groupsE] + int$int
  yhat
  # yhat <- b1[groupsG] + b2[groupsE] + aux[,3]

  #yhat <- b1[groupsG] + b2[groupsE] + model_run_BAMMIT$BUGSoutput$mean$blin
  return(yhat)
}

# indexTest <- list(gtest = data$datTest$df$gen,
#                   etest = data$datTest$df$env)

predyTrain <- predFunc(model = model_run_BAMMIT$BUGSoutput,
                  indexTrain = list(groupsG = groupsG,
                                    groupsE = groupsE),
                  indexTest = list(gtest = data$datTrain$df$gen,
                                   etest = data$datTrain$df$env))
round(caret::RMSE(predyTrain, data$datTrain$df$y), 4)


predyTest <- predFunc(model = model_run_BAMMIT$BUGSoutput,
                       indexTrain = list(groupsG = groupsG,
                                         groupsE = groupsE),
                       indexTest = list(gtest = data$datTest$df$gen,
                                        etest = data$datTest$df$env))
round(caret::RMSE(predyTest, data$datTest$df$y), 4)


round(caret::RMSE(yhat, data$datTest$df$y), 4)
round(caret::RMSE(yhat, data$datTrain$df$y), 4)


# Run AMBARTI -------------------------------------------------------------

library(AMBARTI)

groupsG <- dat$df$gen
groupsE <- dat$df$env

df = dat$df
df$g = as.factor(groupsG)
df$e = as.factor(groupsE)
df$y = df$y
df = df[,-which(colnames(df) %in% c('gen','env','blin'))]
df = as.data.frame(df)
y = df$y
#y <- rnorm(120,10,2)
x <-  df[,-which(colnames(df) == 'y')]
fit.ambarti.ammi = ambarti(x = x, y = y, ntrees = 10, nburn = 50, npost = 100, nsteps = 1)
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


   for (i in 1:N_test) {
    # Model for phenotype
    y_test[i] ~ dnorm(mu[i], sigma^-2)
    mu_test[i] = g[gen_test[i]] + e[env_test[i]] + blin_test[i] #gen and env are parameters
    blin_test[i] = sum(lambda[1:Q] * gamma[gen_test[i],1:Q] * delta[env_test[i],1:Q])

    # Clustering model
    gen_test[i] ~ dcat(piG[i, 1:G1])
    env_test[i] ~ dcat(piE[i, 1:G2])

    # Continuous variables
    g_group_test[i, 1:NG] ~ dmnorm(mu_gen[gen_test[i], 1:NG], sgroup_inv_gen)
    e_group_test[i, 1:NE] ~ dmnorm(mu_env[env_test[i], 1: NE], sgroup_inv_env)

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
model_data <- list(N = N, y = dat$df$y, G1 = usedC, G2 = usedC, NG = NG, NE = NE,
                   sG = diag(NG), sE = diag(NE), Q = dat$Q,
                  env = dat$df$env) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
str(model_data)
# Choose the parameters to watch
model_parameters <- c("g", "e", "gen", "env", "piG", "piE", "mu_env", "mu_gen",
                      "mu", "sigma", "sgroup_gen", 'sgroup_env', 'blin', 'gamma', 'delta', 'lambda')

# Run the model
model_run <- R2jags::jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)


plot(model_run)

round(caret::RMSE(model_run$BUGSoutput$median$mu, dat$df$y), 4)
round(caret::RMSE(model_run$BUGSoutput$median$mu, data$datTest$df$y), 4)


# RUN THE MODELS ----------------------------------------------------------


# Create the data sets

rm(list = ls())

N <- 400
NG <- 2
NE <- 2
sg <- 1
se <- 1
#sgroup <-  diag(2)
sigma <- 1
lambda <- 20

trueC <- 3
G1 <- trueC
G2 <- trueC

set.seed(02)
muG <- matrix(runif(NG*G1, -50,50), ncol = NG, nrow = G1)
muE <- matrix(runif(NE*G2, -50,50), ncol = NE, nrow = G2)
dat <- genData(N = N, G1 = G1, G2 = G2, NG = NG, NE = NE,
               muG = muG, muE = muE, sg = sg, se = se,
               sigma = sigma, lambda = lambda)

train <- sample(1:N, size = 200)
test <- (1:N)[-train]
usedC <-  3

# Run CBAMMIT

# Set up the data
model_data_CBAMMIT <- list(N = length(train), y = dat$df$y[train], G1 = usedC, G2 = usedC,
                           NG = NG, NE = NE, sG = diag(NG), sE = diag(NE), Q = dat$Q,
                           g_group = dat$g_group[train,], e_group = dat$e_group[train,],
                           N_test = length(test), g_group_test = dat$g_group[test,],
                           e_group_test = dat$e_group[test,]) #, alphaG = rep(1,G1), alphaE = rep(1,G2))
str(model_data_CBAMMIT)
# Choose the parameters to watch
model_parameters <- c("g", "e", "gen", "env", "piG", "piE", "mu_env", "mu_gen",
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




# C+BAMMIT ----------------------------------------------------------------

model_data_groupsG <- list(N = length(train), G1 = usedC, NG = NG,
                           g_group = dat$g_group[train,], sG = diag(NG),
                           N_test = length(test), g_group_test = dat$g_group[test,])

model_data_groupsE <- list(N = length(train), G2 = usedC, NE = NE,
                           e_group = dat$e_group[train,], sE = diag(NE),
                           N_test = length(test), e_group_test = dat$e_group[test,])

# Choose the parameters to watch
model_parameters_groupsG <- c("gen", "piG","mu_gen", "sgroup_gen", 'gen_tes')
model_parameters_groupsE <- c("env", "piE", "mu_env", 'sgroup_env', 'env_test')

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

groupsG_test <- as.integer(model_run_groupsG$BUGSoutput$median$gen)
groupsE_test <- as.integer(model_run_groupsE$BUGSoutput$median$env_test)


# Set up the data
model_data_BAMMIT <- list(N = length(train), y = dat$df$y[train],
                          G1 = usedC, G2 = usedC,
                          gen = groupsG, env = groupsE,
                          gen_test = groupsG_test, env_test = groupsE_test,
                          N_test = length(test), Q = dat$Q)

str(model_data_BAMMIT)

# Choose the parameters to watch
model_parameters_BAMMIT <- c("g", "e", "sigma", 'mu', 'mu_test' ,'blin', 'lambda', 'gamma', 'delta')

# Run the model
model_run_BAMMIT <- R2jags::jags(
  data = model_data_BAMMIT,
  parameters.to.save = model_parameters_BAMMIT,
  model.file = textConnection(model_code_BAMMIT)
)

caret::RMSE(dat$df$y[train], model_run_BAMMIT$BUGSoutput$median$mu)
caret::RMSE(dat$df$y[test], model_run_BAMMIT$BUGSoutput$median$mu_test)


df <- data.frame(model = c(rep('CBAMMIT', 18), rep('C+BAMMIT', 18)),
                 data = c(rep('Train', 9), rep('Test', 9), rep('Train', 9), rep('Test', 9)),
                 true = as.factor(c(rep(3,3), rep(10,3), rep(20,3), rep(3,3), rep(10,3), rep(20,3),
                          rep(3,3), rep(10,3), rep(20,3), rep(3,3), rep(10,3), rep(20,3))),
                 used = c(3,10,20, 3,10,20, 3,10,20, 3,10,20,3,10,20,3,10,20,
                          3,10,20,3,10,20,3,10,20,3,10,20,3,10,20,3,10,20),
                 rmse = c(1.1762, 0.4854,0.1746,1.1699,0.9816,0.5256,1.4342,1.3143,1.0872,
                          8.0502, 5.0290, 1.9648, 2.4870, 1.3466, 1.2468, 2.1902, 1.9301,1.7448,
                          1.2949, 0.9912,0.9714, 2.3987,1.3579, 1.3094, 1.9261, 1.7481, 1.5775,
                          9.9169, 9.0173, 8.2034, 3.0735,2.9330, 2.7146, 2.3235, 2.2624, 1.9054))
library(ggplot2)

df |> ggplot(aes(x = used, y = rmse, colour = true)) +
  geom_line()+
  geom_point()+
  scale_color_manual(values=c('steelblue', 'firebrick', 'darkgreen')) +
  facet_grid(data~model) +
  theme_bw()
