library(readxl)
library(dplyr)
library(patchwork)
library(R2jags)
library(readr)
library(mclust)
library(tidyr)

rm(list = ls())

# ireland <- read_excel("/Volumes/ALESSA HD/PhD/bammit/Real data/ireland.xlsx", sheet = "Optimal")
# vcu <- readRDS(("/Volumes/ALESSA HD/PhD/bammit/Real data/vcus/complete_vcu.rds"))
# saveRDS(ireland, file = 'data/ireland.rds')
# saveRDS(vcu, file = 'data/vcu.rds')


# irl <- ireland |> select("Genotype", "Location", "Bloc", 'Year',
#                   "yld_ton_ha", 'plot_length', 'plot_width', 'ave_moist')
# names(irl) <- c("Genotype", "Environment", 'Bloc', 'Year', 'Yield',
#                 'plot_length', 'plot_width', 'ave_moist')
#
# irl <- irl |> mutate_at(c("Genotype", "Environment", 'Bloc'), as.factor)
#
# levels(irl$Genotype) <- paste0("g", 1:length(levels(irl$Genotype)))
# levels(irl$Environment) <- paste0("e", 1:length(levels(irl$Environment)))
# levels(irl$Bloc) <- paste0("b", 1:length(levels(irl$Bloc)))
# head(irl.df, 20)
# irl.df <- as.data.frame(irl)
# names(ireland)
# str(ireland)

ireland <- readRDS("~/Documents/GitHub/clusteringModel/data/ireland.rds")
ireland <-  ireland |>
  filter(Year %in% c(2010)) |>
  dplyr::select("Genotype", "Location", "Bloc", "yld_ton_ha", 'plot_length', 'ave_moist') |>
  rename(Yield = yld_ton_ha, Environment = Location) |>
  mutate_at(c("Yield",'plot_length', 'ave_moist'), as.numeric) |>
  mutate_at(c("Genotype", "Environment", 'Bloc'), as.factor) |>
  group_by(Genotype, Environment) |>
  mutate(yield = mean(Yield),
         length = mean(plot_length),
         moist = mean(ave_moist)) |>
  summarise(yield = sort(unique(yield)),
            length = sort(unique(length)),
            moist = sort(unique(moist))) |>
  ungroup() |>
  dplyr::select("Genotype", "Environment", "yield", 'moist')

# Rename factor levels  to be anonymous
levels(ireland$Genotype) <- paste0("g", 1:length(levels(ireland$Genotype)))
levels(ireland$Environment) <- paste0("e", 1:length(levels(ireland$Environment)))

str(ireland)
hist(ireland$moist)




X <- ireland[,-c(1,2,3,4)]
groups <- ireland$Environment
clPairs(X, groups)
BIC <- mclustBIC(X)
mod1 <- Mclust(X, x = BIC)
par(mfrow = c(1,2))
plot(BIC)
plot(mod1, what = "classification")


# innovar data

innovar_all_data <- read_delim("data/innovar-all-data.csv",
                               delim = ";", escape_double = FALSE, trim_ws = TRUE)
names(innovar_all_data)
str(innovar_all_data)


innovar_all_data |> filter(c(crop == 'Bread Wheat', year == '2020'))

innovar_sel <- unique(innovar_all_data[,c('Yield','variety', 'WD')]) # add: PP, TGW, WD

dat <- innovar_sel |>
  fill(WD, .direction = "up") |>
  drop_na(c(Yield, WD)) |>
  mutate_at(c("variety"), as.factor)

levels(dat$variety) <- paste0("g", 1:length(levels(dat$variety)))
dat$variety
nrow(dat)

dat0 <- dat[dat$WD == 0,]
dat3 <- dat[dat$WD == 3,]
dat9 <- dat[dat$WD == 9,]

dat <- rbind(dat3[1:200,], dat9[1:200,])

dat$Yield |> hist()
dat$variety |> unique()
dat$WD |> hist()

# innovar <- innovar %>%
#   mutate_all(as.numeric) %>%
#   fill(names(.), .direction = "up")
# innovar$PP <- as.factor(innovar$PP)
# innovar$PP[is.na(innovar$PP)] <- '9'
# str(innovar)
# tail(innovar)
# hist(innovar$TGW)


# X <- innovar[,-3]
# groups <- innovar$PP
# clPairs(X, groups)
# BIC <- mclustBIC(X)
# mod1 <- Mclust(X, x = BIC)
# par(mfrow = c(1,2))
# plot(BIC)
# plot(mod1, what = "classification")

# apply model

model_code <- '
model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] =    e[env[i]] +  g[gen[i]] # Note that env here is a parameter but gen is not
    #blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

    # Continuous environmental variable
    t[i] ~ dnorm(mu_env[env[i]], st^-2)
    # Clustering model
    env[i] ~ dcat(pi[1:G])
    #env[i] ~ dcat(pi[i, 1:G])

    # for (j in 1:G) {
    #   exp_theta[i, j] <- exp(theta[i, j])
    #   pi[i, j] <- exp(theta[i, j]) / sum(exp_theta[i, 1:G])
    #   theta[i, j] ~ dnorm(0, 6^-2)
    # }

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
N <- nrow(dat)
y <- dat$Yield
G <- 2
I <- length(levels(dat$variety))
gen <- dat$variety
t <- dat$WD

# innovar data

# N <- nrow(ireland)
# y <- ireland$yield
# G <- 2
# I <- length(levels(ireland$Genotype))
# gen <- ireland$Genotype
# t <- ireland$moist

# vcu data
# N <- nrow(vcu)
# y <- vcu$yield
# G <- 2
# t <- vcu$moist_cont


model_data <- list(N = N, y = y, G = G, I = I, gen = gen,# Q = 1,
                   t = t, alpha = rep(1,G))

# Choose the parameters to watch
model_parameters <- c("g", "e", "env", "pi", "mu_env", "mu", "sigma", "st")
#model_parameters <- c("e", "env", "pi", "mu_env", "mu", "sigma", "st")

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
p1 <- qplot(model_run$BUGSoutput$median$env, ireland$Environment) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_bw() +
  labs(x = 'estimated env', y = 'ireland env')

# Overall predictions
p2 <- qplot(model_run$BUGSoutput$median$mu, y) +
  geom_abline() + theme_bw() +
  labs(x = expression(mu))

p1+p2




# With interaction --------------------------------------------------------


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
N <- nrow(ireland)
y <- ireland$yield
G <- 2
I <- length(levels(ireland$Genotype))
gen <- ireland$Genotype
t <- ireland$moist
Q <- 1

model_data <- list(N = N, y = y, G = G, I = I, gen = gen, Q = Q,
                   t = t, alpha = rep(1,G))

# Parameters to watch
model_parameters <- c("g", "e", "env", "pi", "mu_env", "mu", "sigma", "st", "blin")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

# Results and output of the real data
plot(model_run)
print(model_run)


# Multivariate ------------------------------------------------------------

rm(list = ls())


ireland <- readRDS("~/Documents/GitHub/clusteringModel/data/ireland.rds")
ireland <-  ireland |>
  filter(Year %in% c(2019)) |>
  select("Genotype", "Location", "Bloc", "yld_ton_ha", 'plot_length', 'ave_moist') |>
  rename(Yield = yld_ton_ha, Environment = Location) |>
  mutate_at(c("Yield",'plot_length', 'ave_moist'), as.numeric) |>
  mutate_at(c("Genotype", "Environment", 'Bloc'), as.factor) |>
  group_by(Genotype, Environment) |>
  mutate(yield = mean(Yield),
         length = mean(plot_length),
         moist = mean(ave_moist)) |>
  summarise(yield = sort(unique(yield)),
            length = sort(unique(length)),
            moist = sort(unique(moist))) |>
  ungroup() |>
  select("Genotype", "Environment", "yield", "length", 'moist')

# Rename factor levels  to be anonymous
levels(ireland$Genotype) <- paste0("g", 1:length(levels(ireland$Genotype)))
levels(ireland$Environment) <- paste0("e", 1:length(levels(ireland$Environment)))

str(ireland)
hist(ireland$moist)
hist(ireland$length)


# Set up the data
N <- nrow(ireland)
y <- ireland$yield
NE <- 2
G <- 2
I <- length(levels(ireland$Genotype))
gen <- ireland$Genotype
t <- cbind(ireland$moist, ireland$length)
Q <- 1

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


plot(model_run)
print(model_run)


# Plot the posterior cluster membership
qplot(model_run$BUGSoutput$median$env, dat$df$env) +
  geom_jitter(width = 0.1, height = 0.1) + theme_light() +
  labs(x = 'environment cluster', y = 'environment cluster estimated')

# Overall predictions
qplot(model_run$BUGSoutput$median$mu, dat$df$y) +
  geom_abline() + theme_light() +
  labs(x = expression(mu), y = 'y')

# Prediction of genotype effects
qplot(model_run$BUGSoutput$median$g, dat$g) +
  geom_abline() + theme_light() +
  labs(x = 'genotype estimated', y = 'genotype')

# Prediction of environment effects
qplot(model_run$BUGSoutput$median$e, dat$e) +
  geom_abline() + theme_light() +
  labs(x = 'environment estimated', y = 'environment')

