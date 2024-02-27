# Load necessary libraries
library(dplyr)

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