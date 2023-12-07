# load daily rainfall data data
df2010 <- read.csv("/Users/alaninglis/Desktop/rainfall/2010.csv", sep = ",")
df2011 <- read.csv("/Users/alaninglis/Desktop/rainfall/2011.csv", sep = ",")
df2012 <- read.csv("/Users/alaninglis/Desktop/rainfall/2012.csv", sep = ",")
df2013 <- read.csv("/Users/alaninglis/Desktop/rainfall/2013.csv", sep = ",")
df2014 <- read.csv("/Users/alaninglis/Desktop/rainfall/2014.csv", sep = ",")
df2015 <- read.csv("/Users/alaninglis/Desktop/rainfall/2015.csv", sep = ",")
df2016 <- read.csv("/Users/alaninglis/Desktop/rainfall/2016.csv", sep = ",")
df2017 <- read.csv("/Users/alaninglis/Desktop/rainfall/2017.csv", sep = ",")
df2018 <- read.csv("/Users/alaninglis/Desktop/rainfall/2018.csv", sep = ",")
df2019 <- read.csv("/Users/alaninglis/Desktop/rainfall/2019.csv", sep = ",")


dfList <- list(
  df2010,
  df2011,
  df2012,
  df2013,
  df2014,
  df2015,
  df2016,
  df2017,
  df2018,
  df2019
)

# remove 1st 2 columns
newDfList <- lapply(dfList, function(x) x[, -c(1,2)])

# Calculate column means for each data frame
meanDfList <- lapply(newDfList, colMeans)

# convert to mm
dfListmm <- lapply(meanDfList, function(df) {
  as.data.frame(lapply(df, function(col) round(col/10, 3)))
})


# rename columns
dfListRename <- lapply(dfListmm, function(df) {
  colnames(df) <- gsub(pattern = "X\\d{4}", replacement = "", x = colnames(df))
  return(df)
})

dfListLong <- lapply(dfListRename, function(df) {
  colnames(df) <- gsub(pattern = "(\\d{2})(\\d{2})", replacement = "\\1_\\2", x = colnames(df))
  return(df)
})

#reshape
library(reshape2)
dfReshaped <- lapply(dfListLong, function(df) {
  melted_df <- melt(df, id.vars = NULL, variable.name = "Date", value.name = "Avg_Rainfall")
  return(melted_df)
})



# ------------------------------------------------------------------------

# monthly average ---------------------------------------------------------



# The average_temps list now contains the average temperature for each month
objects_to_remove <- ls()[ls() != "dfReshaped"]
rm(list=objects_to_remove)
rm(objects_to_remove)


# Function to compute monthly average
compute_monthly_avg <- function(df) {
  df %>%
    mutate(monthName = substr(Date, 1, 2)) %>%
    group_by(monthName) %>%
    summarise(avgMonthRain = mean(Avg_Rainfall, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(monthName = month.abb[as.numeric(monthName)])
}

# Apply the function to each data frame in the list
avgMonthRain <- lapply(dfReshaped, compute_monthly_avg)


years <- as.character(2010:2019)
avgMonthRain <- setNames(avgMonthRain, years)

save.image(file = 'avgMonthRain.Rdata')


# -------------------------------------------------------------------------


# name lists
years <- as.character(2010:2019)
dfListFinal <- setNames(dfReshaped, years)

yearlyAvgDf <- lapply(dfListFinal, function(df) {
  avg_rainfall <- mean(df$Avg_Rainfall) 
})

yearlyAvgDf <- unlist(yearlyAvgDf)

avgYearRain <- data.frame(Year = years, avg_rainfall_mm = yearlyAvgDf)

rownames(avgYearRain) <- NULL


save.image(file = 'avgRainfall_mm.Rdata')

# -------------------------------------------------------------------------










# remove garbage
rm(df2010,
   df2011,
   df2012,
   df2013,
   df2014,
   df2015,
   df2016,
   df2017,
   df2018,
   df2019,
   dfList,
   newDfList,
   meanDfList,
   dfListmm,
   dfListRename,
   dfListLong,
   #dfReshaped,
   yearlyAvgDf,
   dfListFinal, 
   merged_avg_rainfall_df,
   years)
