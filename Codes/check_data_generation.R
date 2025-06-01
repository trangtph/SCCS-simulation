# This script is to check if the distribution of the exposure and event in the dataset is as specified


test_data_generation <- function(n = 100000, obs_time = 500, baseline_odds = 2e-5, IRR = 2, risk_length=28, rep)
{
  beta0 <- log(baseline_odds)   
  beta1 <- log(IRR)
  
  dataset <- data.frame()
  
  # Assign exposure status
  exposed <- rbinom(n, 1, 0.8)
  # If exposed, assign exposure date
  exposure_start <- rep(NA, n)
  exposure_start[exposed == 1] <- sample(1:(obs_time - risk_length), 
                                         sum(exposed), replace = TRUE)
  exposure_end <- exposure_start + risk_length - 1
  
  # Expand the data to long format: one row per individual per day
  id <- rep(1:n, each = obs_time)
  day <- rep(1:obs_time, times = n)
  
  # Daily exposure status
  exposure_start_long <- rep(exposure_start, each = obs_time)
  exposure_end_long <- rep(exposure_end, each = obs_time)
  exposure_status <- ifelse(!is.na(exposure_start_long) & day >= exposure_start_long & day <= exposure_end_long, 1, 0)
  # Daily outcome status
  logit_p <- beta0 + beta1 * exposure_status
  p_event <- plogis(logit_p)
  outcome <- rbinom(length(p_event), 1, p_event)
  
  # Long data format
  long_data <- data.frame(
    id = id,
    day = day,
    exposure = exposure_status,
    outcome = outcome,
    exposure_start = exposure_start_long,
    exposure_end = exposure_end_long
  )
  
  
  # Calculate odds of events
  no_exposure_days <- long_data %>% filter(exposure == 0)
  odds_event_no_exposure <- sum(no_exposure_days$outcome)/(nrow(no_exposure_days) - sum(no_exposure_days$outcome))
  exposure_days <- long_data %>% filter(exposure == 1)
  odds_event_exposure <- sum(exposure_days$outcome)/(nrow(exposure_days) - sum(exposure_days$outcome))
  
  #proportion of vaccination
  prop_exposed <- mean(!is.na(long_data$exposure_start[!duplicated(long_data$id)]))
  
  check_data <- data.frame(
    rep = rep,
    odds_event_no_exposure = odds_event_no_exposure,
    odds_event_exposure = odds_event_exposure,
    odds_ratio = odds_event_exposure/ odds_event_exposure,
    prop_exposed = prop_exposed
  )
  
  return(check_data)
}

check_data_output <- test_data_generation(rep=1)

set.seed(222)
nrep = 5
check_data_output2 <- foreach( i = 1:nrep, .combine="rbind") %do% {
  checkdata <- test_data_generation(rep=i)
}
