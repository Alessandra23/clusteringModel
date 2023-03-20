library(R2jags)

rm(list = ls())

genData <- function(N, G1, G2, muG, muE, sg, se, sgroup, sigma){

  g <- rnorm(G1, 0, sg)
  e <- rnorm(G2, 0, se)

  # clusters
  gen <- sample(1:G1, size = N, replace = TRUE)
  env <- sample(1:G2, size = N, replace = TRUE)

  # generate mean of the model
  mu_y <- g[gen] + e[env]
  y <- rnorm(N, mu_y, sigma)
  g_group <- rnorm(N, muG[gen], sgroup)
  e_group <- rnorm(N, muE[env], sgroup)
  df <- data.frame(y = y, env = env, gen = gen, g_group = g_group, e_group = e_group)


  return(list(df = df,
              g = g,
              e = e))

}

G1 <- 3
G2 <- 3
N <- 200
muG <- c(-5, 10, 30)
muE <- c(-10, 0, 10)
sg <- 1
se <- 1
sgroup <-  1
sigma <- 1
set.seed(02)
dat <- genData(N = N, G1 = G1, G2 = G2, muG = muG, muE = muE, sg = sg, se = se, sgroup = sgroup, sigma = sigma)

# hist(dat$df$g_group, breaks = 30, freq = FALSE)
# for (g in 1:G1) curve(dnorm(x, mean = muG[g])/G1, col = g, add = TRUE)
# hist(dat$df$e_group, breaks = 30, freq = FALSE)
# for (g in 1:G2) curve(dnorm(x, mean = muE[g])/G2, col = g, add = TRUE)
# hist(dat$df$y)

#plot env clustering
ggplot(dat$df, aes(x = e_group)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  stat_function(fun = function(x) dnorm(x, mean = muE[1])/G2, color = "steelblue",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muE[2])/G2, color = "firebrick",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muE[3])/G2, color = "chartreuse4",size = 0.7) +
  theme_bw() +
  xlab('Temperature')


#plot gen clustering
ggplot(dat$df, aes(x = g_group)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  stat_function(fun = function(x) dnorm(x, mean = muG[1])/G1, color = "steelblue",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muG[2])/G1, color = "firebrick",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muG[3])/G1, color = "chartreuse4",size = 0.7) +
  theme_bw()

#plot y
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
    mu[i] = g[gen[i]] + e[env[i]]

    # Clustering model
    gen[i] ~ dcat(piG[1:G1])
    env[i] ~ dcat(piE[1:G2])

    # Continuous genotopyc variable
    g_group[i] ~ dnorm(mu_gen[gen[i]], sgroup^-2)

    # Continuous environmental variable
    e_group[i] ~ dnorm(mu_env[env[i]], sgroup^-2)
   }

  # Priors

  # Prior on clusters membership
  piG ~ ddirch(alphaG)
  piE ~ ddirch(alphaE)

  for (g in 1:G1) {
    mu_gen_raw[g] ~ dnorm(0, 100^-2)
  }
  # Make sure these are in order to avoid label switching
  mu_gen <- sort(mu_gen_raw[1:G1])


  for (g in 1:G2) {
    mu_env_raw[g] ~ dnorm(0, 100^-2)
  }
  # Make sure these are in order to avoid label switching
  mu_env <- sort(mu_env_raw[1:G2])

  # Prior on genotype effect
  for(i in 1:G1) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:G2) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sgroup ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)

  }
'

# Set up the data
model_data <- list(N = N, y = dat$df$y, G1 = G1, G2 = G2, gen = dat$df$gen, env = dat$df$enf,
                   g_group = dat$df$g_group, e_group = dat$df$e_group, alphaG = rep(1,G1),
                   alphaE = rep(1,G2))

# Choose the parameters to watch
model_parameters <- c("g", "e", "gen", "env", "piG", "piE", "mu_env", "mu_gen", "mu", "sigma", "sgroup")

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
qplot(model_run$BUGSoutput$median$gen, dat$df$gen) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_light() +
  xlab('Estimated group') +
  ylab('True group')

qplot(model_run$BUGSoutput$median$env, dat$df$env) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_light() +
  xlab('Estimated group') +
  ylab('True group')


# Overall predictions
qplot(model_run$BUGSoutput$median$mu, dat$df$y) + geom_abline() + theme_light()

# Prediction of genotype effects
qplot(model_run$BUGSoutput$median$g, dat$g) + geom_abline() + theme_light()

# Prediction of environment effects
qplot(model_run$BUGSoutput$median$e, dat$e) + geom_abline() + theme_light()

