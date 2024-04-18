### Aquatic Forecast Workflow ###
#devtools::install_github("eco4cast/neon4cast")
library(tidyverse)
#library(neon4cast)
library(lubridate)
#install.packages("rMR")
#library(rMR)

forecast_date <- Sys.Date()
noaa_date <- Sys.Date() - days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet

#Step 0: Define team name and team members 
team_info <- list(team_name = "AquaticEcosystems",
               team_list = list(
                  list(individualName = list(givenName = "Sebastian", 
                                        surName = "Sanchez"),
                      organizationName = "Boston University",
                      electronicMailAddress = "bastian6@bu.edu"),
                  list(individualName = list(givenName = "Trenton", 
                                        surName = "Mulick"),
                      organizationName = "Boston University",
                      electronicMailAddress = "tmulick@bu.edu"),
                  list(individualName = list(givenName = "Jianrui", 
                                        surName = "Liu"),
                      organizationName = "Boston University",
                      electronicMailAddress = "jianrui@bu.edu"),
                  list(individualName = list(givenName = "Matthew", 
                                        surName = "Manberg"),
                      organizationName = "Boston University",
                      electronicMailAddress = "mmanber1@bu.edu"),
                  list(individualName = list(givenName = "Alejandro", 
                                        surName = "Diaz"),
                      organizationName = "Boston University",
                      electronicMailAddress = "alejdiaz@bu.edu"))
)

## Load required functions
if(file.exists("site_temp_oxygen_data.R"))      source("site_temp_oxygen_data.R")
# if(file.exists("NOAA_Forecast_download.R"))    source("NOAA_Forecast_download.R")
if(file.exists("DLM.R"))    source("DLM.R")

# DLM model

dlm <- DLM_function()

### Step 1: Download Required Data
download <- download_targets()       ## Y variables

target <- download[[1]]
site_names <- download[[2]]

site_data  <- download_site_meta()
ox_temp <- filter_ox_temp(target)

target_oxygen <- ox_temp[[1]]
target_temp <- ox_temp[[2]]

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
