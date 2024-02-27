### Aquatic Forecast Workflow ###
devtools::install_github("eco4cast/neon4cast")
library(tidyverse)
library(neon4cast)
library(lubridate)
install.packages("rMR")
library(rMR)

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
                      electronicMailAddress = "lejdiaz@bu.edu"))
)

## Load required functions
if(file.exists("site_temp_oxygen_data.R"))      source("site_temp_oxygen_data.R")


### Step 1: Download Required Data
target     <- download_targets()       ## Y variables
site_data  <- download_site_meta()
target     <- merge_met_past(target)   ## append met data (X) into target file
met_future <- download_met_forecast(forecast_date) ## Weather forecast (future X)

# ## visual check of data
# ggplot(target, aes(x = temperature, y = air_temperature)) +
#   geom_point() +
#   labs(x = "NEON water temperature (C)", y = "NOAA air temperature (C)") +
#   facet_wrap(~site_id)

# met_future %>% 
#   ggplot(aes(x = datetime, y = air_temperature, group = parameter)) +
#   geom_line() +
#   facet_grid(~site_id, scale ="free")
