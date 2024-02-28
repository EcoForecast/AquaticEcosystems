# Load necessary libraries
library(dplyr)
install.packages("neonUtilities")
library(neonUtilities)

#' Download target data
#'
#' This function downloads the target data from a specific URL and removes any rows with missing values.
#' It also extracts the unique site IDs from the target data.
#'
#' @return A list containing the target data and the unique site IDs.
download_targets <- function(){
  target <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)
  target <- na.omit(target) # Get days with no observations
  site_names <- unique(target$site_id)
  return(list(target, site_names))
}
# What is this URL and how is it different from what is at the bottom?


#' Download site metadata
#'
#' This function downloads the site metadata from a specific URL and filters it to include only aquatic sites.
#'
#' @return The filtered site metadata.
download_site_meta <- function(){
  site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
    filter(aquatics == 1)
  return(site_data)
}

#' Filter target data for oxygen and temperature
#'
#' This function filters the target data to include only rows where the variable is "oxygen" or "temperature".
#'
#' @param target The target data.
#' @return A list containing the filtered target data for oxygen and temperature.
filter_ox_temp <- function(target){
  target_oxygen <- target |> 
    filter(variable == "oxygen")
  target_temp <- target |> 
    filter(variable == "temperature")
  return(list(target_oxygen, target_temp))
}





site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  filter(aquatics == 1)
sites <- site_data$field_site_id

# Download past year of historical data for all sites via API request
# Data includes water quality (chl a and dissolved oxygen), temp at specific depth of lakes, and temp at surface of streams
# We may need to get an API token to help with this, keep getting rate limited
# Each API request requires the user to input a "y" or "n" in the console to confirm the download, issue for automating?

# Water quality is largest data product by a lot, >1GB, split into two requests, only Jan - Apr
water_quality_2023 <- loadByProduct(dpID="DP1.20288.001", site=sites, package="basic", startdate = "2023-01", enddate = "2023-04", tabl = "waq_instantaneous")
#water_quality_2024 <- loadByProduct(dpID="DP1.20288.001", site=sites, package="basic", startdate = "2024-01", enddate = "2024-02", tabl = "waq_instantaneous", include.provisional = T)
 # Above line not working, skipping it for now

temp_at_depth <- loadByProduct(dpID="DP1.20264.001", site=sites, package = "basic", startdate = "2023-01", enddate = "2024-02", tabl = "TSD_30_min")

temp_at_surface <- loadByProduct(dpID="DP1.20053.001", site=sites, package = "basic", startdate = "2023-01", enddate = "2024-02", tabl = "TSW_30min")


# Plotting one site for illustration instead of all 34 for each variable
# Temp at Depth
temp_ad_tbl <- temp_at_depth$TSD_30_min |> 
  select(siteID, startDateTime, endDateTime, tsdWaterTempMean)

temp_ad_tbl <- na.omit(temp_ad_tbl)

temp_ad_site1 <- temp_ad_tbl |> 
  filter(siteID == sites[2])

plot(temp_ad_site1$startDateTime, temp_ad_site1$tsdWaterTempMean,
     pch = 16,
     cex = 0.5,
     main = sites[2],
     ylab = "Temperature At Depth (C)",
     xlab = "Date"
     )

# Water Quality, Dissolved Oxygen & Chlorophyll a
water_quality_site1 <- water_quality$waq_instantaneous |> 
  filter(siteID == sites[2]) |> 
  select(siteID, startDateTime, dissolvedOxygen, chlorophyll)

water_quality_site1 <- na.omit(water_quality_site1)

plot(water_quality_site1$startDateTime, water_quality_site1$dissolvedOxygen,
     ylab = "Dissolved Oxygen",
     xlab = "Date",
     main = sites[2],
     type = 'l'
     )
plot(water_quality_site1$startDateTime, water_quality_site1$chlorophyll,
     ylab = "Chlorophyll",
     xlab = "Date",
     main = sites[2],
     type = 'l'
     )

# Temp at Surface
temp_as_site1 <- temp_at_surface$TSW_30min |> 
  filter(siteID == sites[1]) |> 
  select(siteID, startDateTime, surfWaterTempMean)

temp_as_site1 <- na.omit(temp_as_site1)

plot(temp_as_site1$startDateTime, temp_as_site1$surfWaterTempMean,
     xlab = "Date",
     ylab = "Temperature at Surface (C)",
     main = sites[1],
     type = 'l'
     )
