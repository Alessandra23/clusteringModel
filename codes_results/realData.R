library(R2jags)

# Real data

# Loading weather data
load("/Users/alessalemos/Documents/GitHub/clusteringModel/data/weatherData/avgTemps.Rdata")
load("/Users/alessalemos/Documents/GitHub/clusteringModel/data/weatherData/avgRainfall_mm.Rdata")

plot(avgTemp2010$AvgTemperature)
avgTemp <- lapply(avgTemps, function(x){
  mean(x[,2])
}) |> tibble::as_tibble() |> reshape2::melt()
colnames(avgTemp) <- c('Year', 'Temperature')

avg_weather <- cbind(avgTemp, Rainfall = avgYearRain$avg_rainfall_mm)

hist(avg_weather$Rainfall)
hist(avg_weather$Temperature)
plot(avg_weather$Rainfall, avg_weather$Temperature)

avg <- avg_weather[-1, ]
# clustering <- Mclust(avg_weather[-1, -1], G = 2)
#avg_weather$Cluster <- as.factor(clustering$classification)

clustering <- Mclust(avg[, -1], G = 2)
avg$Cluster <- as.factor(clustering$classification)

avg |>
  ggplot(aes(Rainfall, Temperature, color = Cluster)) +
  geom_point(size = 2.5) +
  theme_bw(base_size = 16) +
  scale_color_manual(values=c("firebrick", "steelblue"))

avg_weather |>
  ggplot(aes(Rainfall, Temperature, colour = Year)) +
  geom_point(size = 2.5) +
  theme_bw(base_size = 16)



# Join with he Irish data

load("~/Documents/GitHub/bammit/Real data/ireland.RData")
ireland_weather <- merge(ireland, avg_weather, by = 'Year')


# Applying the model


cbammit_real_data <- function(data, data_test, Q = 1, mmu = 10, smu = 2, stheta = 1, gamma = 10000,
                              a = 0.1, b = 0.1, nthin = 1, nburnin = 2000){

  Q <- Q
  Y <- data$Yield
  N <- length(data$Yield)
  B1 <- length(levels(unique(data$Genotype)))
  B2 <- length(levels(unique(data$Environment)))
  # B3 <- length(levels(unique(data$Year)))
  B4 <- length(levels(unique(data$Bloc)))

  var1 <- data$Genotype
  var2 <- data$Environment
  # var3 <- data$Year
  var4 <- data$Bloc

  Y_test <- data_test$Yield
  N_test <- length(data_test$Yield)
  var1_test <- data_test$Genotype
  var2_test <- data_test$Environment
  # var3_test <- data_test$Year
  var4_test <- data_test$Bloc

  model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[n] =  b1[var1[n]] + b2[var2[n]] + b3[var3[n]] + b4[var2[n], var4[n]] + int[n] # clustering just var3 (year)
     int[n] = sum(lambda[1:Q]*beta1[var1[n],1:Q]*beta2[var2[n],1:Q]*beta3[var3[n], 1:Q])

    # # Clustering model
    t[i, 1:NE] ~ dmnorm(muT[env[i], 1:NE], st_inv)
    var3[i] ~ dcat(pi[i, 1:G])

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



}

