library(readr)
library(dplyr)
library(lubridate)
library(ecoforecastR)  # Ensure this package is installed and loaded

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
  model <- list(obs = "y", fixed = "~ 1 + X", n.iter = 5000)

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
##' @param  M   = model matrix
##' @param  mu0 = initial condition mean vector
##' @param  P0  = initial condition covariance matrix
##' @param  Q   = process error covariance matrix
##' @param  R   = observation error covariance matrix
##' @param  Y   = observation matrix (with missing values as NAs), time as col's
##'
##' @return list
##'  mu.f, mu.a  = state mean vector for (a)nalysis and (f)orecast steps
##'  P.f, P.a    = state covariance matrix for a and f
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
    print(t)
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

##' Kalman Filter: Analysis step
##' @param  mu.f = Forecast mean (vector)
##' @param  P.f  = Forecast covariance (matrix)
##' @param  Y    = observations, with missing values as NAs) (vector)
##' @param  R    = observation error covariance (matrix)
##' @param  H    = observation matrix (maps observations to states)
KalmanAnalysis <- function(mu.f, P.f, Y, R, H, I) {
    # Identify which observations are not missing
    obs = !is.na(Y)
    
    # Check the length of obs against R to handle dimensionality correctly
    if (length(obs) != ncol(R)) {
        stop("Mismatch in dimension of R and expected observations")
    }
    
    if (any(obs)) {
        print("Observations found, updating estimates...")
        # Adjust H and R for observed values
        H_obs = H[obs, , drop = FALSE]
        R_obs = R[obs, obs, drop = FALSE]  # Correctly subset R based on observed data

        # Kalman Gain computation
        S = H_obs %*% P.f %*% t(H_obs) + R_obs
        K = P.f %*% t(H_obs) %*% solve(S)

        # Update estimates
        mu.a = mu.f + K %*% (Y[obs] - H_obs %*% mu.f)
        P.a = (I - K %*% H_obs) %*% P.f
    } else {
        # If no observations, updates are not possible
        print("No observations found, maintaining previous estimates...")
        mu.a = mu.f
        P.a = P.f
    }

    return(list(mu.a = mu.a, P.a = P.a))
}


##' Kalman Filter: Forecast Step
##' @param mu.a = analysis posterior mean (vector)
##' @param P.a  = analysis posterior covariance (matrix)
##' @param M    = model (matrix)
##' @param  Q   = process error covariance (matrix)
KalmanForecast <- function(mu.a, P.a, M, Q) {
  # Ensure mu.a is a column matrix
  if (is.vector(mu.a)) {
    mu.a <- matrix(mu.a, nrow = length(mu.a), ncol = 1)
  }

  # Check dimensions to prevent non-conformable errors
  if (!all(dim(M) == c(nrow(mu.a), nrow(mu.a)))) {
    stop("Dimension mismatch: M does not match the dimensions of mu.a")
  }
  if (!all(dim(P.a) == c(nrow(mu.a), nrow(mu.a)))) {
    stop("Dimension mismatch: P.a dimensions are incorrect")
  }
  if (!all(dim(Q) == dim(M))) {
    stop("Dimension mismatch: Q does not match the dimensions of M")
  }

  # Forecast mean and covariance computation
  mu.f = M %*% mu.a
  P.f = Q + M %*% P.a %*% t(M)

  return(list(mu.f = mu.f, P.f = P.f))
}




# Assuming the DLM_function and Kalman filter functions are defined as provided previously

# Run the DLM function to get initial parameters
results <- DLM_function()

# Print the summary of DLM parameters
print(results$params)

# Assuming Y is a vector of observations (nrow x 1)
Y = matrix(results$Y, nrow = length(results$Y), ncol = 1)
M = diag(1, length(results$mu0))  # State transition matrix
Q = diag(0.01, length(results$mu0))  # Process noise
R = diag(1, length(results$Y))  # Observation noise matrix, adjust based on actual observation count

kf_results = KalmanFilter(M, results$mu0, results$P0, Q, R, Y)

# Extracting and printing the results for analysis
print(kf_results$mu.a)  # Analysis mean states
print(kf_results$P.a)   # Analysis covariance matrices


