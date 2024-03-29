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
#devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
library(ecoforecastR)
library(zoo)
library(padr)

target <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)
oxygen <- target |> 
  dplyr::filter(variable == "oxygen")


site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv")


site_data <- site_data |> 
  dplyr::filter(aquatics == 2)

site_ids <- site_data$field_site_id

# Load historical data using neon4cast
df_past <- neon4cast::noaa_stage3()


# Jags Random Walk Model
RandomWalk = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"

# Filter data for only one site
site_ox <- oxygen |> 
  dplyr::filter(site_id == "ARIK")

y <- site_ox$observation
data <- list(y=y,n=length(y),      ## data
             x_ic=log(1000),tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
             )

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(y.samp)),  ## initial guess on process precision
                    tau_obs=5/var(y.samp))        ## initial guess on obs precision
}

j.model   <- jags.model (file = textConnection(RandomWalk),
                             data = data,
                             inits = init,
                             n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                                n.iter = 1000)
plot(jags.out)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                                n.iter = 10000)

out <- as.matrix(jags.out)         ## convert from coda to matrix  

time <- site_ox$datetime
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x

ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],
     type='n',
     ylim=range(y,na.rm=TRUE),
     ylab="Oxygen"
     )
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)