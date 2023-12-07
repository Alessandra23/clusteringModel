# #### reading the genomic kernel... basically a "distance matrix" between genotypes. 
# data_gbkernel <- read.table("HEL/GENOTYPE-GBkernel.txt", header = TRUE)
# data_gkkernel <- read.table("HEL/GENOTYPE-GBkernel.txt", header = TRUE)
# 
# # ## This dataset comes from https://data.cimmyt.org/file.xhtml?fileId=5604&version=1.2
# # Read the file into a data frame
# ## year of experiments = 2015. 
# data <- read.table("HEL/pheno_gy_geno.txt", header = TRUE)
# # If your file is comma-separated, you can use read.csv() instead:
# # data <- read.csv(file_path)
# # Print the first few rows of the data frame
# head(data)
# unique(data$Location)## contains 5 environments 
# ## we need to discover from the paper which are the environments:
# # Load the ggplot2 package
# # Load the ggplot2 package
# library(ggplot2)
# 
# # Assuming your data frame is named "df"
# ggplot(data, aes(x = factor(Location), y = Y)) +
#   geom_boxplot() +
#   labs(x = "Location", y = "Y") +
#   ggtitle("Boxplots of Y Across Locations")
# 
# 
# ## From crossing graphs of page 1997 in 
# ## SOUZA 2017 we can see that the levels of 
# ## envs are : e1=IP, e2:NM, e3=PM, e4=SE,
# ## e5=SO
# 
# # lets double check with EH and PH
# 
# data_2 <- read.table("HEL/pheno_geno_ph_eh.txt", header = TRUE)
# # If your file is comma-separated, you can use read.csv() instead:
# # data <- read.csv(file_path)
# # Print the first few rows of the data frame
# head(data_2)
# unique(data_2$Location)## contains 5 environments 
# ## we need to discover from the paper which are the environments:
# # Load the ggplot2 package
# # Load the ggplot2 package
# library(ggplot2)
# ## Y1 is plant height  
# # Assuming your data frame is named "df"
# ggplot(data_2, aes(x = factor(Location), y = Y1)) +
#   geom_boxplot() +
#   labs(x = "Location", y = "Y") +
#   ggtitle("Boxplots of Y Across Locations")
# 
# # Y2 is ear height 
# ggplot(data_2, aes(x = factor(Location), y = Y2)) +
#   geom_boxplot() +
#   labs(x = "Location", y = "Y") +
#   ggtitle("Boxplots of Y Across Locations")
# 
# data$Location <- factor(data$Location)  # Convert Location to a factor variable
# data_2$Location <- factor(data_2$Location)  # Convert Location to a factor variable
# 
# # Define a mapping for Location codes to actual location names
# location_mapping <- c("Ipiaçu", "NovaMutum", "Pato de Minas", "Sertanópolis", "Sorriso")
# location_mapping_2 <- c("Ipiaçu","Pato de Minas", "Sertanópolis")
# 
# # Update the Location column with the actual location names
# data$Location_name <- location_mapping[data$Location]
# data_2$Location_name <- location_mapping_2[data_2$Location]
# 
# ###
# colnames(data)[3]="Yield"
# colnames(data_2)[c(3,4)]=c("plant_height", "ear_height")
# head(data)
# head(data_2)
# 
# head(data %>% arrange(Germplasm_id))
# head(data_2 %>% arrange(Germplasm_id))
# # Merge the two datasets based on Germplasm_id and Location_name
# # Merge the datasets by Location, Germplasm_id, and Location_name
# merged_data <- merge(data, data_2, by = c("Location", "Germplasm_id", "Location_name"), all = TRUE)
# head(merged_data %>% arrange(Germplasm_id))
# 
# write_rds(merged_data,"HEL/phenomic-data.rds")
# ## Now with the information about the enviroments we can adapt the dataset
# # ##  
# # e1= Ipiaçu, 
# # e2=Nova Mutum,
# # e3=Patos de Minas, 
# # e4=Sertanópolis,
# # e5=Sorriso
# # Assuming your data frame is named "data"
# 
# 
# # # Define the location names and corresponding environment vectors
# # location_names <- c("Ipiaçu", "Nova Mutum", "Patos de Minas", "Sertanópolis", "Sorriso")
# # latitude_vector <- c(18.419, 13.059, 18.349, 23.039, 12.329)
# # longitude_vector <- c(-49.569, -56.059, -46.319, -51.029, -55.429)
# # height_vector <- c(452, 460, 832, 361, 365)
# # 
# # # Create the correspondence data frame
# # correspondence_data <- data.frame(
# #   Location = location_names,
# #   Lat = latitude_vector,
# #   Long = longitude_vector,
# #   Height_above_sea = height_vector
# # )
# 
# 
# # Load the required packages
# library(rnaturalearth)
# library(ggplot2)
# library(viridis)
# library(sf)
# 
# # Define the data frame with location information
# correspondence_data <- data.frame(
#   Location = c("Ipiaçu", "Nova Mutum", "Patos de Minas", "Sertanópolis", "Sorriso"),
#   Lat = c(18.419, 13.059, 18.349, 23.039, 12.329),
#   Long = c(-49.569, -56.059, -46.319, -51.029, -55.429),
#   Height_above_sea = c(452, 460, 832, 361, 365)
# )
# 
# # Get a world map with the boundaries of Brazil using rnaturalearth
# world <- ne_countries(country = "Brazil", returnclass = "sf")
# 
# # Create a base map of Brazil
# brazil_map <- ggplot(data = world) +
#   geom_sf(fill = "lightgray", color = "gray") +
#   coord_sf(xlim = c(-75, -30), ylim = c(-35, 10), expand = FALSE) +
#   labs(title = "Brazil Map") +
#   theme_void()
# 
# # Add the cities to the map
# brazil_map +
#   geom_point(data = correspondence_data, aes(x = Long, y = Lat, color = Height_above_sea), size = 3) +
#   geom_text(data = correspondence_data, aes(x = Long, y = Lat, label = Location), size = 3, color = "black") +
#   scale_color_viridis(name = "Height Above Sea Level") +
#   theme_void()
# 
# 
# 
# # Install and load the nasapower package
# 
# #install.packages("nasapower", dependencies = TRUE)
# #install.packages("nasapower")
# # Define the location (latitude and longitude)
# # Specify the variable of interest (e.g., temperature)
# 
# # # Download the data
# # ag_d <- get_power(
# #   community = "ag",
# #   lonlat = c(151.81, -27.48),
# #   pars = c("RH2M", "T2M", "PRECTOTCORR"),
# #   dates = "1985-01-01",
# #   temporal_api = "daily"
# # )
# # 
# # ag_d
# 
# # Define the variable of interest (e.g., temperature)
# variables <- c(
#   "T2M_MAX", "T2M_MIN", "T2M", "T2M_RANGE", "PRECTOTCORR",
#   "RH2M", "WS10M", "WD10M",  "ALLSKY_SFC_SW_DWN"
# )
# 
# 
# # Initialize an empty list to store data for each location
# data_list <- list()
# 
# # Loop through each location and download data
# for (i in 1:nrow(correspondence_data)) {
#   location <- correspondence_data[i, ]
#   data_wth <- nasapower::get_power(
#     community = "ag",  # Use "ag" community for agriculture
#     lonlat = c(location$Long, location$Lat),
#     pars = variables,
#     dates = c("2014-09-01","2016-04-31"),  # Specify the date range as needed
#     temporal_api = "monthly"  # Use "daily" for daily data
#   )
#   data_list[[location$Location]] <- data_wth
# }
# 
# # Print the first few rows of data for each location
# for (i in 1:length(data_list)) {
#   cat("Location:", names(data_list)[i], "\n")
#   print(head(data_list[[i]]))
# }

#write_rds(data_list, "HEL/weather_data.rds")




# # View the resulting dataset
# unique(data$Germplasm_id)
# 
# 
# library(tidyr)
# 
# # Create an empty list to store individual data frames
# df_list <- list()
# 
# # Loop through the list and extract the relevant data
# for (location_name in names(data_list)) {
#   location_data <- data_list[[location_name]]
#   
#   # Gather the data from wide to long format
#   location_df <- gather(location_data, Month, Value, -c(LON, LAT, PARAMETER, YEAR))
#   
#   # Extract desired parameters
#   parameters <- c(
#     "T2M", "RH2M", "WD10M", "WS10M", "T2M_MAX",
#     "T2M_MIN", "T2M_RANGE", "PRECTOTCORR"
#   )
#   
#   # Filter the dataframe to include only desired parameters
#   location_df <- location_df[location_df$PARAMETER %in% parameters, ]
#   
#   # Rename columns
#   colnames(location_df) <- c("LON", "LAT", "PARAMETER", "YEAR", "Month", parameters)
#   
#   # Extract year and location name
#   year <- location_df$YEAR
#   location <- rep(location_name, length(year))
#   
#   # Add Year and Location columns
#   location_df$Year <- year
#   location_df$Location <- location
#   
#   # Append the dataframe to the list
#   df_list[[location_name]] <- location_df
# }
# 
# # Combine all dataframes into a single dataframe
# final_df <- do.call(rbind, df_list)
# 
# # Print the first few rows of the final dataframe
# head(final_df)
# View(final_df)
# 
# library(tidyr)
# 
# # List of parameters to pivot
# parameters <- c(
#   "T2M", "RH2M", "WD10M", "WS10M", "T2M_MAX",
#   "T2M_MIN", "T2M_RANGE", "PRECTOTCORR"
# )
# 
# # Initialize an empty list to store the pivoted dataframes
# pivoted_data_list <- list()
# library(tidyverse)
# # Loop through each parameter and pivot the dataset
# for (param in parameters) {
#   pivoted_data <- final_df %>%
#     filter(PARAMETER == param) %>%
#     dplyr::select(-PARAMETER) %>%
#     spread(param, T2M)
#   
#   # Rename columns as needed
#   colnames(pivoted_data)[1:3] <- c("lon", "lat", "location")
#   
#   # Add the pivoted dataframe to the list
#   pivoted_data_list[[param]] <- pivoted_data
# }
# 
# # Access each parameter's pivoted dataframe like pivoted_data_list$T2M, pivoted_data_list$RH2M, etc.
