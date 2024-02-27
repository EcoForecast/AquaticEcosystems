# This script downloads and creates scatter plots for the historical 
# temperature and oxygen data for all Aquatic Ecosystems NEON field sites

library(dplyr)

target <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)

site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  filter(aquatics == 1)

target <- na.omit(target) # Get days with no observations

# Filter for just oxygen measurements
target_oxygen <- target |> 
  filter(variable == "oxygen")

# Filter for just temperature measurements
target_temp <- target |> 
  filter(variable == "temperature")

site_names <- unique(target$site_id)

# Plot oxygen and temperature for each site
for(site in site_names){
  dates_oxygen <- target_oxygen |> 
    filter(site_id == site) |> 
    select(datetime)
  site_oxygen <- target_oxygen |> 
    filter(site_id == site) |> 
    select(observation)
  
  dates_temp <- target_temp |> 
    filter(site_id == site) |> 
    select(datetime)
  site_temp <- target_temp |> 
    filter(site_id == site) |> 
    select(observation)
  
  plot(dates_oxygen$datetime, site_oxygen$observation,
       xlab = "Date",
       ylab = "Oxygen",
       main = site,
  )
  plot(dates_temp$datetime, site_temp$observation,
       xlab = "Date",
       ylab = "Temperature",
       main = site,
  )
}