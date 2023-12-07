# Code for paper

# Real data


# https://figshare.com/articles/dataset/Additional_file_2_of_Bayesian_modelling_of_phosphorus_content_in_wheat_grain_using_hyperspectral_reflectance_data/22601378


# devtools::install_github('allogamous/EnvRtype',force=TRUE)
# library(envRtype)
library(agridat)

agridat

load("/Users/alessalemos/Downloads/maizeYield.RData")

maizeYield |> names()

data('wheat', package = 'agridat')

dat <- agridat::ars.earlywhitecorn96
plot(dat$flower, dat$moisture)
plot(dat$flower)
plot(dat$moisture)

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
Gused <- 2
model_data <- list(N = N, y = dat$df$y, G = Gused, I = I, NE = NE, gen = dat$df$gen,
                   t = dat$t, stE = diag(NE), Q = dat$Q)#, alpha = rep(1,G))

# Choose the parameters to watch
model_parameters <- c("g", "e", "env", "pi", "muT", "mu", "sigma", "st", 'blin', 'gamma',
                      'delta', 'lambda')
#model_parameters <- c("g")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)


plot(model_run)



## real data

weather_data <- readRDS("~/Documents/GitHub/clusteringModel/data/danilo_alessa_final/weather_data.rds")
phenomic_data <- readRDS("~/Documents/GitHub/clusteringModel/data/danilo_alessa_final/phenomic-data.rds")



weather_data$Ipiaçu |> names()


weather_data$Ipiaçu$ANN |> hist(xlab = 'variable', main = ' ')

