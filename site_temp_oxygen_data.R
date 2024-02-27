# This script downloads and creates scatter plots for the historical 
# temperature and oxygen data for all Aquatic Ecosystems NEON field sites

library(dplyr)



download_targets <- function(){
  target <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)
  target <- na.omit(target) # Get days with no observations
  site_names <- unique(target$site_id)
  return(list(target, site_names))
}

download_site_meta <- function(){
  site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
    filter(aquatics == 1)
  return(site_data)
}


filter_ox_temp <- function(target){
  target_oxygen <- target |> 
    filter(variable == "oxygen")
  target_temp <- target |> 
    filter(variable == "temperature")
  return(list(target_oxygen, target_temp))
}

