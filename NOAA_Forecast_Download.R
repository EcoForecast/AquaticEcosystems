library(neon4cast)
library(dplyr)
library(lubridate)
library(ggplot2)
library(arrow)
library(glue)
library(readr)

target <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)

site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  filter(aquatics == 1)

site_ids <- target |> na.omit() |> distinct(site_id, variable) |> 
  filter(variable %in% c("oxygen", "temperature")) |>
  count(site_id) |> filter(n==2) |> pull(site_id)

# Load historical data using neon4cast
df_past <- neon4cast::noaa_stage3()

# Helper function for historical data
noaa_mean_historical <- function(df_past, site, var) {
  df_past |>
    dplyr::filter(site_id == site, variable == var) |>
    dplyr::rename(ensemble = parameter) |>
    dplyr::select(datetime, prediction, ensemble) |>
    dplyr::mutate(date = as_date(datetime)) |>
    dplyr::group_by(date) |>
    dplyr::summarize(mean_prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(datetime = date) |>
    dplyr::mutate(mean_prediction = if_else(var == "air_temperature", mean_prediction - 273.15, mean_prediction)) |>
    dplyr::collect()
}

# Helper function to download future forecast data for a site and variable
noaa_mean_forecast <- function(site, var, reference_date) {
  endpoint <- "data.ecoforecast.org"
  bucket <- glue::glue("neon4cast-drivers/noaa/gefs-v12/stage2/parquet/0/{reference_date}")
  s3 <- arrow::s3_bucket(bucket, endpoint_override = endpoint, anonymous = TRUE)
  
  ds <- arrow::open_dataset(s3)
  
  forecast_date <- as_datetime(reference_date)
  
  # Filter, select, and process data
  ds %>%
    dplyr::filter(site_id == site,
                  datetime >= forecast_date,
                  variable == var) %>%
    dplyr::select(datetime, prediction, parameter) %>%
    dplyr::mutate(datetime = as_date(datetime)) %>%
    dplyr::group_by(datetime, parameter) %>%
    dplyr::summarize(mean_prediction = mean(prediction, na.rm = TRUE), .groups = "drop") %>%
    dplyr::select(datetime, mean_prediction, parameter) %>%
    dplyr::rename(ensemble = parameter) %>%
    dplyr::collect()
}

################################################################################
#Historical Forecast

# Variables of interest
variable_of_interest <- "air_temperature" # Adjust as needed

# Initialize lists to store data frames for "air_temperature"
past_data_list <- list()

# Loop through site IDs directly
for (site_id in site_ids) {
  # Load historical data using neon4cast
  df_past <- neon4cast::noaa_stage3()
  
  # Fetch and process historical data for the site
  df_past_processed <- noaa_mean_historical(df_past, site_id, variable_of_interest)
  past_data_list[[site_id]] <- df_past_processed
}

# Save the past data for "air_temperature" across all sites
save(past_data_list, file = "past_data_air_temperature.Rdata")

################################################################################
#Future Forecast

variable_of_interest <- c("air_temperature") # Adjust as needed
future_data_list <- list()
for (site_id in site_ids) {
  # Forecast date for future data
  forecast_date <- Sys.Date() - lubridate::days(50) # Adjust as needed, the latest date not sure
  
  # Process future forecast data for the site
  df_future <- noaa_mean_forecast(site_id, variable_of_interest, forecast_date)
  future_data_list[[site_id]] <- df_future
}
save(future_data_list, file = paste0("future_data_", variable_of_interest, ".Rdata"))
