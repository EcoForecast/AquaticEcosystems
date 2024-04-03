library(neon4cast)
library(dplyr)
library(lubridate)
library(ggplot2)
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

site = "BARC"

site_ox <- oxygen |> 
  dplyr::filter(site_id == site)
site_temp <- temperature |> 
  dplyr::filter(site_id == site)

y <- site_ox$observation
temp <- site_temp$observation

# Assuming site_ox and site_temp are data frames with a datetime column
merged_data <- merge(site_ox, site_temp, by = "datetime", all = FALSE)

time <- merged_data$datetime

# Now, y and temp are vectors of the same length
y <- merged_data$observation.x
temp <- merged_data$observation.y

data <- list(y=y,n=length(y),      ## data
             x_ic=log(1000),tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1,           ## process error prior
             temp = temp
)

#model = list(obs="y",fixed="~ 1 + X + temp", n.iternumber = 20000)

model = list(obs="y",fixed="~ 1 + X", n.iternumber = 20000)

ef.out <- ecoforecastR::fit_dlm(model=model,data)

params <- window(ef.out$params,start=1000)
plot(params)
summary(params)
cor(as.matrix(params))
pairs(as.matrix(params))

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


forecastN <- function(IC ,betaIntercept, betaX ,Q=0,n=Nmc){

  NT = 30
  N <- matrix(NA,n,NT)  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:NT){
    mu = Nprev + betaX * Nprev + betaIntercept ## mean
    N[,t] <- rnorm(n,mu,Q)                         ## predict next step
    Nprev <- N[,t]                                  ## update IC
  }
  return(N)
}

## parameters
params <- as.matrix(ef.out$params)
param.mean <- apply(params,2,mean)

IC <- as.matrix(ef.out$predict)


forecast <- forecastN(mean(IC[,"x[1651]"]), param.mean[1], param.mean[2], Q = 0, n = 1)
NT = 30
time1 = 1:NT


x = 1:length(ci[2,])

time2 = length(x):(length(x)+NT-1)

plot(x,ci[2,],
     type='n',
     ylim=range(y,na.rm=TRUE),
     ylab="Dissolved Oxygen",
     xlim = c(1500,length(ci[2,]) + NT),
     )
ecoforecastR::ciEnvelope(x,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(x,y,pch="+",cex=0.5)

lines(time2,forecast,col="purple",lwd=3)