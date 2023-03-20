library(nasaweather)
library(nasadata)
library(nasapower)

atmos <- nasaweather::atmos
nasaweather

library(devtools)
install_github("nasa/NASAaccess", build_vignettes = TRUE)
library(NASAaccess)


nasaweather::atmos

dat <- atmos[, c('temp', 'pressure')]


atmos$lat == 53.3

hist(atmos$lat)

DFworkbook2020
DFworkbook2021

atmos |> filter(str_detect(lat, "^53"))

dat$temp |> hist()
dat$pressure |> hist()


hist(DFworkbook2020)



