

# Initialize an empty list to store the reshaped data frames for each month
list_data <- list()


# Loop through each of the 12 files
for (month in 1:12) {
  # Construct the file name based on the month using the provided naming convention
  file_name <- paste0("~/Desktop/temperature/min/2019/IRL_DLY_TN_2019", sprintf("%02d", month), "_grid_IDW_09Z.csv")
  
  # Read the data set
  data <- read.csv(file_name)
  
  # Drop the first two columns
  data <- data[, -(1:2)]
 
  
  # Calculate the average temperature for the month
  avg_temp <- colMeans(data, na.rm = TRUE)
  
  # Reshape the data so that each day is a separate row
  long_data <- reshape2::melt(avg_temp)
  
  # Create a data frame with the average temperature and the temperature values
  month_data <- data.frame(date = rownames(long_data), temperature = long_data$value)
  
  # Append to the list
  list_data[[month]] <- month_data
}

# -------------------------------------------------------------------------


# Assuming your list of data frames is named 'df_list'

#library(lubridate)
#library(dplyr)

# Function to process a single data frame
process_df <- function(df) {
  # Extract year, month, and day
  df$Year <- as.numeric(substr(df[,1], 2, 5))
  df$Month <- as.numeric(substr(df[,1], 6, 7))
  df$Day <- as.numeric(substr(df[,1], 8, 9))
  
  # Convert to Date class
  df$Date <- as.Date(with(df, paste(Year, Month, Day, sep="-")), format="%Y-%m-%d")
  
  # Get month name
  df$MonthName <- months(df$Date, abbreviate = FALSE)
  
  # Group by month name and calculate average temperature
  df_avg <- df %>%
    group_by(MonthName) %>%
    summarise(AvgTemperature = mean(df[,2], na.rm = TRUE))
  
  return(df_avg)
}

# Apply the function to each data frame in the list
new_df_list <- lapply(list_data, process_df)

dfListmm <- lapply(new_df_list, function(df) {
  as.data.frame(df)
})



temp2019MonthMin <- setNames(dfListmm, month.abb)

# -------------------------------------------------------------------------

objects_to_remove <- ls()[ls() != "temp2019MonthMin"]
rm(list=objects_to_remove)
rm(objects_to_remove)



save.image(file = 'temp2019MonthMin.Rdata')








