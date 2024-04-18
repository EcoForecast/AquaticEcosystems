library(readr)
library(dplyr)
library(lubridate)
library(ecoforecastR)  # Ensure this package is installed and loaded
library(ggplot2)

DLM_function <- function() {
  # Load target data
  target <- read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)
  site_data <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv")

  # Filter for aquatic sites and specific variables
  site_data <- site_data |> filter(aquatics == 1)
  site = "BARC"
  oxygen <- target |> filter(variable == "oxygen", site_id == site)

  # Assuming oxygen data is sorted or needs sorting by datetime for sequential analysis
  oxygen <- oxygen |> arrange(datetime)
  y <- oxygen$observation
  time <- oxygen$datetime

  # Prepare data for DLM
  data <- list(y = y, n = length(y))

  # Define the DLM model; simplified model since `X` isn't used
  model <- list(obs = "y", fixed = "~ 1 + X", n.iter = 2000)

  # Run the DLM
  ef.out <- ecoforecastR::fit_dlm(model = model, data = data)

  # Extract and process results from DLM
  params <- window(ef.out$params, start = 1000)
  summary_params <- summary(params)

  # Initial conditions for Kalman Filter
  mu0 <- as.numeric(head(y, 1))  # First observation of y
  P0 <- diag(10, length(mu0))    # Large initial covariance

  # Return relevant data and DLM results
  return(list(
    Y = y,                       # Observational data for Kalman filter if needed
    mu0 = mu0,                   # Initial mean for Kalman filter
    P0 = P0,                     # Initial covariance matrix for Kalman filter
    time = time,                 # Time vector
    DLM_results = ef.out,        # DLM fit results
    params = summary_params      # Summary of parameters from DLM
  ))
}


##'  Kalman Filter

KalmanFilter <- function(M,mu0,P0,Q,R,Y){
  
  ## storage
  nstates = nrow(Y)  
  nt = ncol(Y)
  mu.f  = matrix(NA,nstates,nt+1)  ## forecast mean for time t
  mu.a  = matrix(NA,nstates,nt)  ## analysis mean for time t
  P.f  = array(NA,c(nstates,nstates,nt+1))  ## forecast variance for time t
  P.a  = array(NA,c(nstates,nstates,nt))  ## analysis variance for time t

  ## initialization
  mu.f[,1] = mu0
  P.f[,,1] = P0
  I = diag(1,nstates)

  ## run updates sequentially for each observation.
  for(t in 1:nt){

    ## Analysis step: combine previous forecast with observed data
    KA <- KalmanAnalysis(mu.f[,t],P.f[,,t],Y[,t],R,H=I,I)
    mu.a[,t] <- KA$mu.a
    P.a[,,t] <- KA$P.a
    
    ## Forecast step: predict to next step from current
    KF <- KalmanForecast(mu.a[,t],P.a[,,t],M,Q)
    mu.f[,t+1] <- KF$mu.f
    P.f[,,t+1] <- KF$P.f
  }
  
  return(list(mu.f=mu.f,mu.a=mu.a,P.f=P.f,P.a=P.a))
}

KalmanAnalysis <- function(mu.f,P.f,Y,R,H,I){
  obs = !is.na(Y) ## which Y's were observed?
  if(any(obs)){
    H <- H[obs,]                                              ## observation matrix
    K <- P.f %*% t(H) %*% solve(H%*%P.f%*%t(H) + R[obs,obs])  ## Kalman gain
    mu.a <- mu.f + K%*%(Y[obs] - H %*% mu.f)                  ## update mean
    P.a <- (I - K %*% H)%*%P.f                                ## update covariance
    ## Note: Here's an alternative form that doesn't use the Kalman gain
    ## it is less efficient due to the larger number of matrix inversions (i.e. solve)
    ## P.a <- solve(t(H)%*%solve(R[obs,obs])%*%(H) + solve(P.f))                             
    ## mu.a <- P.a %*% (t(H)%*%solve(R[obs,obs])%*%Y[obs] + solve(P.f)%*%mu.f)
  } else {
    ##if there's no data, the posterior is the prior
    mu.a = mu.f
    P.a = P.f
  }
  return(list(mu.a=mu.a,P.a=P.a))
}


KalmanForecast <- function(mu.a,P.a,M,Q){
  mu.f = M%*%mu.a
  P.f  = Q + M%*%P.a%*%t(M)
  return(list(mu.f=mu.f,P.f=P.f))
}




# Run the DLM function to get initial parameters
results <- DLM_function()

# Define the number of states based on your data
nstates <- length(results$Y)

# Log transformation of data
Y <- results$Y
Y <- matrix(Y, ncol = 1)  # Ensure Y is a column matrix for consistency in Kalman Filter

# Define a simple state transition matrix M, assuming simple evolution without spatial interactions
alpha <- 0.05
M <- diag(1 - alpha, nstates) + alpha * diag(1, nstates)

# Define the error covariance matrices with correct dimensions
tau_proc <- rep(0.01, nstates)
Q <- diag(tau_proc)
tau_obs <- 0.1
R <- diag(tau_obs, nstates)

# Initial conditions based on historical data or estimated from DLM
mu0 <- results$mu0
P0 <- results$P0
if (is.null(dim(P0))) {
    P0 <- diag(0.01, nstates)
}

# Run the Kalman Filter
KF_results <- KalmanFilter(M, mu0, P0, Q, R, Y)

# Plot the actual and predicted oxygen levels
time <- results$time
actual <- results$Y
predicted <- KF_results$mu.a

df <- data.frame(Time = time, Actual = actual, Predicted = predicted)
ggplot(df, aes(x = Time)) +
  geom_line(aes(y = Actual), color = "blue") +
  geom_line(aes(y = Predicted), color = "red") +
  labs(y = "Oxygen Level", 
       title = "Actual vs Predicted Oxygen Levels",
       subtitle = "Blue: Actual, Red: Predicted")