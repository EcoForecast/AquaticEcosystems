library(neon4cast)
library(dplyr)
library(lubridate)
library(arrow)
library(glue)
library(readr)
library(rjags)
#library(rnoaa)
library(daymetr)
# devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
library(ecoforecastR)
library(zoo)
library(padr)


DLM_function <- function(site) {

target <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)
oxygen <- target |> 
  dplyr::filter(variable == "oxygen")

temperature <- target |> 
  dplyr::filter(variable == "temperature")

site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv")


site_data <- site_data |> 
  dplyr::filter(aquatics == 1)

site_ids <- site_data$field_site_id

# Load historical data using neon4cast
# df_past <- neon4cast::noaa_stage3()

#site = "ARIK"

site_ox <- oxygen |> 
  dplyr::filter(site_id == site)
site_temp <- temperature |> 
  dplyr::filter(site_id == site)

y <- site_ox$observation
temp <- site_temp$observation


### Weather Forecast Data ###

### Weather Hindcast Data ###
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

df_past <- neon4cast::noaa_stage3()

# Fetch and process historical data for the site
df_past_processed <- noaa_mean_historical(df_past, site, "air_temperature")
# Merge in past NOAA data into the targets file, matching by date.
site_target <- target |>
  dplyr::select(datetime, site_id, variable, observation) |>
  dplyr::filter(variable %in% c("temperature", "oxygen"), 
                site_id == site_id) |>
  tidyr::pivot_wider(names_from = "variable", values_from = "observation") |>
  dplyr::left_join(df_past_processed, by = c("datetime"))


# Assuming site_ox and site_temp are data frames with a datetime column
merged_data <- merge(site_ox, df_past_processed, by = "datetime", all = FALSE)

time <- merged_data$datetime

# Now, y and temp are vectors of the same length
y <- merged_data$observation
temp <- merged_data$mean_prediction

data <- list(y = y,n = length(y),
             x_ic=mean(y, na.rm = T),
             tau_ic= 1/sd(y),
             a_obs=1,
             r_obs=1,          
             a_add=1,
             r_add=1,           
             temp = temp
)

#model = list(obs="y",fixed="~ 1 + X + temp", n.iternumber = 20000)

model = list(obs="y",fixed="~ 1 + X + temp", n.iternumber = 10000)

ef.out <- ecoforecastR::fit_dlm(model=model,data)

params <- window(ef.out$params,start=1000)
plot(params)
summary(params)
#cor(as.matrix(params))
#pairs(as.matrix(params))

out <- as.matrix(ef.out$predict)
ci <- apply(out,2,quantile,c(0.025,0.5,0.975))
plot(time,ci[2,],
     type='n',
     ylim=range(y,na.rm=TRUE),
     ylab="Dissolved Oxygen",
     )
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)



#### Milestone 6 ####


forecastN <- function(IC, betaIntercept, betaX, betaTemp, temp , Q = 0, n = Nmc){

  NT = 30
  N <- matrix(NA,n,NT)  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:NT){
    mu = Nprev + betaX * Nprev + betaTemp * temp[t] + betaIntercept  ## mean
    N[,t] <- rnorm(n,mu,Q)                         ## predict next step
    Nprev <- N[,t]                                  ## update IC
  }
  return(N)
}

## parameters

reference_date <- Sys.Date() - 1

df_future <- neon4cast::noaa_stage2(start_date = as.character(reference_date))

noaa_air_temp_future <- df_future |> 
  dplyr::filter(datetime >= lubridate::as_datetime(reference_date), 
                variable %in% c("air_temperature")) |>
  dplyr::filter(site_id == site) |> 
  dplyr::rename(ensemble = parameter) |> 
  dplyr::select(datetime, prediction, ensemble) |> 
  dplyr::mutate(date = as_date(datetime)) |> 
  dplyr::group_by(date) |> 
  dplyr::summarize(mean_prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  dplyr::rename(datetime = date) |> 
  dplyr::collect() |> 
  dplyr::mutate(mean_prediction = mean_prediction - 273.15)


params <- as.matrix(ef.out$params)
param.mean <- apply(params,2,mean)

IC <- as.matrix(ef.out$predict)


#forecast <- forecastN(mean(IC[,"x[1035]"]), param.mean[1], param.mean[2], param.mean[3], noaa_air_temp_future, Q = 0, n = 1)
NT = 30
time1 = 1:NT
x = 1:length(ci[2,])

time2 = length(x):(length(x)+NT-1)

Nmc = 1000
prow = sample.int(nrow(params),Nmc,replace=TRUE)
forecast <- forecastN(IC[prow, ncol(IC)], params[prow,"betaIntercept"], params[prow,"betaX"], params[prow, "betatemp"], noaa_air_temp_future$mean_prediction, Q = 0, n = Nmc)

par(mfrow=c(1,1))

plot(x,ci[2,],
     type='n',
     ylim=range(y,na.rm=TRUE),
     ylab="Dissolved Oxygen",
     # xlim = c(1500,length(ci[2,]) + NT),
     xlim = c(min(time2) - 30, max(time2)),
     main = "Forecast with IC and Parameter Uncertainty"
)
ecoforecastR::ciEnvelope(x,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(x,y,pch="+",cex=0.5)

N.IP.ci = apply(forecast,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha("lightBlue",1))
lines(time2,N.IP.ci[2,],lwd=0.5)
legend("bottomright", legend = c("Data", "CI", "IC Uncertainty", "Parameter Uncertainty"), lty = c(NA,1,1,1), col = c("black", "lightblue", "black", "red"), pch = c("+", NA, NA, NA), cex = 0.7)


###### Organize Output into EFI Standard ######

forecast_date <- Sys.Date()

forecast_dates <- seq(Sys.Date()-1, by = "day", length.out = 30)
forecast_efi <- data.frame(prediction = as.vector(N.IP.ci[2,]))
#cbind(forecast_efi,forecast)
forecast_efi <- forecast_efi |> 
  dplyr::mutate(reference_datetime = forecast_date,
                datetime = forecast_dates,
                site_id = site,
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





##########

#Qmc <- 1/sqrt(params[prow,"Q"])


### Uncertainty Analysis
### calculation of variances
#varI     <- apply(N.I,2,var)
#varIP    <- apply(N.IP,2,var)
#varMat   <- rbind(varI,varIP)

## in-sample stacked area plot
#V.pred.rel.in <- apply(varMat[-5,],2,function(x) {x/max(x)})
#plot(time2,V.pred.rel.in[1,],ylim=c(0,1),type='n',main="Relative Variance: In-Sample",ylab="Proportion of Variance",xlab="time")

#N.cols <- c("red","blue")

#ciEnvelope(time2,rep(0,ncol(V.pred.rel.in)),V.pred.rel.in[1,],col=N.cols[1])
#ciEnvelope(time2,V.pred.rel.in[1,],V.pred.rel.in[2,],col=N.cols[2])
#legend("topleft",legend=c("Process","InitCond"),col=rev(N.cols[-5]),lty=1,lwd=5)
}