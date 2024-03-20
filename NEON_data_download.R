# Load necessary libraries
install.packages("package_name", repos = "https://cloud.r-project.org/")
library(dplyr)
install.packages("neonUtilities")
library(neonUtilities)

# Function to download site data
download_site_data <- function() {
    site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
        filter(aquatics == 1)
    sites <- site_data$field_site_id
    return(list(site_data, sites))
}

# Function to download water quality data
download_water_quality <- function(sites) {
    water_quality_2023 <- loadByProduct(dpID="DP1.20288.001", site=sites, package="basic", startdate = "2023-01", enddate = "2023-04", tabl = "waq_instantaneous")
    return(water_quality_2023)
}

# Function to download temperature at depth data
download_temp_at_depth <- function(sites) {
    temp_at_depth <- loadByProduct(dpID="DP1.20264.001", site=sites, package = "basic", startdate = "2023-01", enddate = "2024-02", tabl = "TSD_30_min")
    return(temp_at_depth)
}

# Function to download temperature at surface data
download_temp_at_surface <- function(sites) {
    temp_at_surface <- loadByProduct(dpID="DP1.20053.001", site=sites, package = "basic", startdate = "2023-01", enddate = "2024-02", tabl = "TSW_30min")
    return(temp_at_surface)
}


# Usage:
site_data_and_sites <- download_site_data()
sites <- site_data_and_sites[[2]]

water_quality_2023 <- download_water_quality(sites)
temp_at_depth <- download_temp_at_depth(sites)
temp_at_surface <- download_temp_at_surface(sites)
