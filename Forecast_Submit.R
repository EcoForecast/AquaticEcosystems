


Submit_function <- function(forecast_df){

###### Organize Output into EFI Standard ######

forecast_date <- Sys.Date()

#forecast_dates <- seq(forecast_date-1, by = "day", length.out = 30)

#forecast_efi <- data.frame(prediction = as.vector(N.IP.ci[2,]))

#cbind(forecast_efi,forecast)

forecast_efi <- forecast_df |> 
  dplyr::mutate(reference_datetime = forecast_date,
                #datetime = forecast_dates,
                #site_id = site,
                variable = "oxygen",
                family = "ensemble",
                duration = "P1D",
                project_id = "neon4cast",
                model_id = "AquaticEcosystemsOxygen",
                parameter = 1,
  ) |> 
  dplyr::select(datetime, reference_datetime, site_id, family, variable, prediction, duration, project_id, model_id, parameter)

team_info <- list(team_name = "AquaticEcosystems",
                  team_list = list(list(individualName = list(givenName = "Trenton", 
                                                              surName = "Mulick"),
                                        organizationName = "Boston University",
                                        electronicMailAddress = "tmulick@bu.edu"),
                                   list(individualName = list(givenName = "Sebastian", 
                                                              surName = "Sanchez"),
                                        organizationName = "Boston University",
                                        electronicMailAddress = "bastian6@bu.edu"),
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
                                        electronicMailAddress = "alejdiaz@bu.edu")
                  )
)

forecast_file <- paste0("aquatics","-",min(forecast_efi$datetime),"-",team_info$team_name,".csv.gz")

write_csv(forecast_efi, forecast_file)
forecast_output_validator(forecast_file)

model_metadata = list(
  forecast = list(
    model_description = list(
      forecast_model_id =  system("git rev-parse HEAD", intern=TRUE), ## current git SHA
      name = "Dynamic Linear Model Oxygen", 
      type = "process",  
      repository = "https://github.com/EcoForecast/AquaticEcosystems"   ## put your REPO here *******************
    ),
    initial_conditions = list(
      status = "absent"
    ),
    drivers = list(
      status = "absent"
    ),
    parameters = list(
      status = "data_driven",
      complexity = 2 # slope and intercept (per site)
    ),
    random_effects = list(
      status = "absent"
    ),
    process_error = list(
      status = "absent"
    ),
    obs_error = list(
      status = "absent"
    )
  )
)


#metadata_file <- neon4cast::generate_metadata(forecast_file, team_info$team_list, model_metadata)

neon4cast::submit(forecast_file = forecast_file, ask = FALSE)


}