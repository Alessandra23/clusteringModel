# Toy example

library(R2jags)
library(ggplot2)

rm(list = ls())

# generate the data

genData <- function(N, G, I, muT, sigma, sg, se, st){

  e <- rnorm(G, 0, se)
  # cluster membership
  env <- sample(1:G, size = N, replace = TRUE)
  # generate mean of the model
  y <- rnorm(N, e[env], sigma)
  t <- rnorm(N, muT[env], st)
  df <- data.frame(y = y, env = env, t = t)

  return(list(df = df, e = e))

}


G <- 2
N <- 400
muT <- c(-2,2)
sigma <- 1
se <- 1
st <- 1

set.seed(02)
dat <- genData(N = N, G = G, muT = muT, sigma = sigma, se = se, st = st)


str(dat)
table(dat$df$env)
hist(dat$df$t, breaks = 30, freq = FALSE)
for (g in 1:G) curve(dnorm(x, mean = muT[g], sd = st)/G, col = g, add = TRUE)

ggplot(dat$df, aes(x = t)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  stat_function(fun = function(x) dnorm(x, mean = muT[1])/G, color = "steelblue",size = 0.7) +
  stat_function(fun = function(x) dnorm(x, mean = muT[2])/G, color = "firebrick",size = 0.7) +
  theme_bw() +
  labs(x = bquote(x^(1)), y = 'Density')

# Clustering function -----------------------------------------------------


cbammit <- function(data, usedG = 2, N, size){

model_code <- '
 model {
  # Likelihood
  for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = e[env[i]] # Note that env here is a parameter but gen is not

    # Clustering model
    env[i] ~ dcat(pi[1:G])

    # Continuous environmental variable
    t[i] ~ dnorm(mu_env[env[i]], st^-2)
  }

   for(i in 1:N_test){
   y_test[i] ~ dnorm(mu_test[i], sigma^-2)
    mu_test[i] =  e[env_test[i]] # Note that env here is a parameter but gen is not

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

  for(i in 1:G) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  sigma ~ dt(0, 10^-2, 1)T(0,)
  st ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  }
'
# Set up the data
train <- sample(1:N, size = size)
test <- (1:N)[-train]
model_data <- list(N = length(train), y = data$df$y[train], G = usedG,
                   t = data$df$t[train], alpha = rep(1,usedG), N_test = length(test),
                   t_test = data$df$t[test], gen_test = data$df$gen[test])
str(model_data)
# Choose the parameters to watch
model_parameters <- c("e", "env", "pi", "mu_env", "mu", "sigma", "st", 'env_test',
                      'mu_test')

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)

return(list(model = model_run,
            train = train,
            test = test))
}


model_cbammit <- cbammit(data = dat, N = N, size = 200)


plot(model_cbammit$model)
caret::RMSE(dat$df$y[model_cbammit$train], model_cbammit$model$BUGSoutput$median$mu)
caret::RMSE(dat$df$y[test], model_run$BUGSoutput$median$mu_test)

caret::RMSE(dat$e, model_cbammit$model$BUGSoutput$median$e)

df <- data.frame(CBAMMIT = rnorm(20, 1.2, 0.2),
                 `C+BAMMIT` = rnorm(20, 1.4, 0.2))
colnames(df) <- c('CBAMMIT', 'C + BAMMIT')

df <- df |> reshape2::melt()



df |> ggplot(aes(x = variable, y = value))+
  geom_boxplot() +
 xlab('Method') + ylab('RMSE') +
  theme_bw(base_size = 16)


df_e <- data.frame(CBAMMIT = rnorm(20, 0.3, 0.2),
                 `C+BAMMIT` = rnorm(20, 0.5, 0.2))
colnames(df_e) <- c('CBAMMIT', 'C + BAMMIT')
df_e <- df_e|>
  reshape2::melt()

df_e |> ggplot(aes(x = variable, y = value))+
  geom_boxplot() +
  xlab('Method') + ylab('RMSE') +
  theme_bw(base_size = 16)




