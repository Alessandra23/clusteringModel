

load("/Users/alaninglis/Desktop/temperature/max/2019/temp2019MonthMax.Rdata")
load("/Users/alaninglis/Desktop/temperature/min/2019/temp2019MonthMin.Rdata")



# Initialize a list to store the average temperatures
avgTemp <- list()

# Loop through each month and compute the average temperature
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
for (month in months) {
  avg_temp <- (temp2019MonthMax[[month]]$AvgTemperature + temp2019MonthMin[[month]]$AvgTemperature) / 2
  avgTemp[[month]] <- data.frame(MonthName = month, AvgTemperature = avg_temp)
}

# Convert the list to a single data frame
avgTemp2019 <- do.call(rbind, avgTemp)

# Reset the row names to avoid repetition
rownames(avgTemp2019) <- NULL


# The average_temps list now contains the average temperature for each month
objects_to_remove <- ls()[ls() != "avgTemp2019"]
rm(list=objects_to_remove)
rm(objects_to_remove)



save.image(file = 'avgTemp2019.Rdata')


# ------------------------------------------------------------------------

# Join all the list objects into a single list


load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2010.Rdata")
load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2011.Rdata")
load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2012.Rdata")
load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2013.Rdata")
load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2014.Rdata")
load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2015.Rdata")
load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2016.Rdata")
load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2017.Rdata")
load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2018.Rdata")
load("/Users/alaninglis/Desktop/temperature/avgTemp/avgTemp2019.Rdata")


# Join all the list objects into a single list
all_temps <- list(avgTemp2010, 
                  avgTemp2011, 
                  avgTemp2012, 
                  avgTemp2013, 
                  avgTemp2014,
                  avgTemp2015, 
                  avgTemp2016, 
                  avgTemp2017, 
                  avgTemp2018, 
                  avgTemp2019)


# name lists
years <- as.character(2010:2019)
avgTemps <- setNames(all_temps, years)



save.image(file = 'avgTemps.Rdata')




