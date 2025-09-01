###################################
### Script for                  ###
### simulating and analysing    ###
### data.                       ###
###################################

#Load packages and data -------------------------------------------------------
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
pacman::p_load(
  foreach,     # foreach loop
  stats,
  extraDistr,
  dplyr,
  SCCS,
  tictoc,      # Measure performance time
  rio,         # Export file  
  here,
  data.table,
  magrittr,     # To use the pipe %>%
  doRNG,        # Reproducible parallel session
  doFuture,
  RhpcBLASctl   # Control threads in parallel session
  
) 

options(scipen = 999)

# Function: generate cohort data -----------------------------------------------

# Four variables: exposure E, outcome Y, confounder C of E and Y that is affected
# by past E, confounder U of C and Y

cohort_time_var_past_e_u <-function(n = 100000, obs_time = 500, risk_window = chosen_risk_window,
                                    p_U = 2e-3, U_on_C = 5, U_on_Y = 5,     
                                    p_C = 1e-3, C_on_C = 0.3, E_on_C = 0.2,
                                    p_E = 3.2e-3, C_on_E = 0.25, E_on_E = 0.001,
                                    baseline_Y = 2e-5, IRR_E = 2, IRR_C = 5)
{
  log_p_U <- log(p_U)
  log_p_C <- log(p_C)
  log_p_E <- log(p_E)
  log_U_on_C <- log(U_on_C)
  log_U_on_Y <- log(U_on_Y)
  log_C_on_C <- log(C_on_C)
  log_E_on_C <- log(E_on_C)
  log_C_on_E <- log(C_on_E)
  log_E_on_E <- log(E_on_E)
  log_baseline_Y <- log(baseline_Y)
  log_IRR_E <- log(IRR_E)
  log_IRR_C <- log(IRR_C)    
  
  U <- base::matrix(0L, nrow = n, ncol = obs_time)
  C <- base::matrix(0L, nrow = n, ncol = obs_time)
  E <- base::matrix(0L, nrow = n, ncol = obs_time)
  Y <- base::matrix(0L, nrow = n, ncol = obs_time)
  
  # Rolling sum of history of U, C and E
  
  roll_U_onC <- base::integer(n)
  roll_U_onY <- base::integer(n)
  roll_C_onC <- base::integer(n)
  roll_E_onC <- base::integer(n)
  roll_C_onE <- base::integer(n)
  roll_E_onE <- base::integer(n)
  roll_C_onY <- base::integer(n)
  roll_E_onY <- base::integer(n)
  
  
  for (t in 1:obs_time) {
    # Remove U, C or E events that fall outside the look-back window 
    # (30 days for C on E, C on Y, /risk window/ for E on Y,
    # 90 days for U on C, U on Y, and 91 days for C on C and E on C)
    if (t > 30) {
      roll_C_onE <- roll_C_onE - C[, t- 30]
      roll_C_onY <- roll_C_onY - C[, t - 30]
      
    }
    if (t > risk_window) {
      roll_E_onY <- roll_E_onY - E[, t - risk_window]
    }
    if (t >90) {
      roll_U_onC <- roll_U_onC - U[, t - 90]
      roll_U_onY <- roll_U_onY - U[, t - 90]       
    }
    if (t > 91) {
      roll_C_onC <- roll_C_onC - C[, t - 91]
      roll_E_onC <- roll_E_onC - E[, t - 91]
    }
    
    # Generate U
    prob_U <- stats::plogis(log_p_U)
    U_t <- stats::rbinom(n, 1, prob_U)
    U[, t] <- U_t
    
    # Update rolling of U
    roll_U_onC <- roll_U_onC + U_t
    roll_U_onY <- roll_U_onY + U_t
    
    # Generate C (depends on U[t-89, t], C[t-90, t-1] and E[t-90, t-1])
    recent_U_onC <- roll_U_onC > 0 # Check recent history of U  
    recent_C_onC <- roll_C_onC > 0 # Check recent history of C
    recent_E_onC <- roll_E_onC > 0 # Check recent history of E
    prob_C <- stats::plogis(log_p_C + log_U_on_C*recent_U_onC + log_C_on_C*recent_C_onC + log_E_on_C*recent_E_onC)
    C_t <- stats::rbinom(n, 1, prob_C)
    C[, t] <- C_t
    
    # Update rolling of C 
    roll_C_onC <- roll_C_onC + C_t
    roll_C_onE <- roll_C_onE + C_t
    roll_C_onY <- roll_C_onY + C_t
    
    # Generate E (depends on C[t-29, t] and E[1, t-1])
    recent_C_onE <- roll_C_onE > 0
    recent_E_onE <- roll_E_onE > 0
    prob_E <- stats::plogis(log_p_E + log_C_on_E*recent_C_onE + log_E_on_E*recent_E_onE)
    E_t <- stats::rbinom(n, 1, prob_E)
    E[, t] <- E_t
    
    # Update rolling of E
    roll_E_onE <- roll_E_onE + E_t
    roll_E_onC <- roll_E_onC + E_t
    roll_E_onY <- roll_E_onY + E_t
    
    
    # Generate Y (depends on U[t-89, t] C[t-29, t], E[t-risk_window+1, t])
    recent_U_onY <- roll_U_onY > 0
    recent_C_onY <- roll_C_onY > 0
    recent_E_onY <- roll_E_onY > 0
    
    prob_Y <- stats::plogis(log_baseline_Y + log_U_on_Y*recent_U_onY + log_IRR_E*recent_E_onY + log_IRR_C*recent_C_onY)
    Y[, t] <- stats::rbinom(n, 1, prob_Y)
    
  }
  
  # Reshape to long format
  
  data_long <- data.table::data.table(
    id = rep(1:n, each = obs_time),
    day = rep(1:obs_time, times = n),
    U = as.vector(t(U)), #transpose vector
    C = as.vector(t(C)),
    E = as.vector(t(E)),
    Y = as.vector(t(Y))
  )
  
  return(data_long)
}

# Function: Reshape data to SCCS-compatible format -------------------------------
# Only keep data on date of occurence of exposure, outcome and covariates

SCCS_reformat_confound <- function(data, risk_window = chosen_risk_window){
  Y_data <- data[Y == 1, .(id, Y_day = day)] #Filter rows with Y==1, select id and day columns
  Y_data[, .row := seq_len(.N), by = id] # For each individual, create event index
  
  E_data <- data[E == 1, .(id, E_day = day)]
  E_data[, .row := seq_len(.N), by = id]
  
  C_data <- data[C == 1, .(id, C_day = day)]
  C_data[, .row := seq_len(.N), by = id]
  
  merge_E_C <- data.table::merge.data.table(E_data, C_data, by = c("id", ".row"), all = TRUE)
  merge_Y_E_C <- data.table::merge.data.table(Y_data, merge_E_C, by = "id", all = TRUE)
  
  # Keep only rows with observed events
  merge_Y_E_C <- merge_Y_E_C[!is.na(Y_day)]
  
  # Add observation window and end of risk periods
  merge_Y_E_C[, `:=`(
    obs_start = 1,
    obs_end = 500,
    E_end = E_day + risk_window - 1,
    C_end = C_day + 30 - 1
  )]
  
  return(as.data.frame(merge_Y_E_C))
}

# Function: analyse data ------------------------------------------------------
# Two approach: ignore C or adjust for C

analyse_sccs_confound <- function(data, rep, C_adjust = TRUE)
{
  if(C_adjust){
    # Fit SCCS model
    model <- SCCS::standardsccs(event ~ E_day + C_day, 
                          indiv= id,         # subject ID
                          astart = obs_start,# start of observation period 
                          aend = obs_end,    # end of observation period 
                          aevent = Y_day,    # event time
                          adrug = cbind(E_day, C_day), # start of exposure
                          aedrug = cbind(E_end, C_end),  # end of risk period
                          expogrp = list(0, 0),       # start of risk period counted from 'adrug'
                          data = data)
    
    # Extract results from SCCS model
    est_E <- stats::coef(model)[1,1]
    se_E <- stats::coef(model)[1,3]
    IRR_E <- exp(est_E)
    IRR_E_CI_Lower <- model$conf.int[1,3]
    IRR_E_CI_Upper <- model$conf.int[1,4]
    est_C <- stats::coef(model)[2,1]
    IRR_C <- stats::coef(model)[2,2]
    IRR_C_CI_Lower <- model$conf.int[2,3]
    IRR_C_CI_Upper <- model$conf.int[2,4]
    n_event <- model$nevent
  }
  else{
    model <- SCCS::standardsccs(event ~ E_day, 
                          indiv= id,      # subject ID
                          astart = obs_start,# start of observation period 
                          aend = obs_end,    # end of observation period 
                          aevent = Y_day, # event time
                          adrug = E_day, # start of exposure
                          aedrug = E_end,  # end of risk period
                          expogrp = 0,       # start of risk period counted from 'adrug'
                          data = data)
    # Extract results from SCCS model
    est_E <- stats::coef(model)[1,1]
    se_E <- stats::coef(model)[1,3]
    IRR_E <- exp(est_E)
    IRR_E_CI_Lower <- model$conf.int[1,3]
    IRR_E_CI_Upper <- model$conf.int[1,4]
    est_C <- NA
    IRR_C <- NA
    IRR_C_CI_Lower <- NA
    IRR_C_CI_Upper <- NA
    n_event <- model$nevent
  }

  res <- data.frame(rep = rep, 
                    est_E, se_E, IRR_E, IRR_E_CI_Lower, IRR_E_CI_Upper, 
                    est_C, IRR_C, IRR_C_CI_Lower, IRR_C_CI_Upper, 
                    n_event, 
                    row.names = NULL 
                    )
  return(res)
}

# Function: bias quantification ------------------------------------------------

bias_quantification <- function(true_IRR_E, true_IRR_C, result_table)
{
  true_gamma1 = log(true_IRR_E)
  true_gamma2 = log(true_IRR_C)
  n_sim = nrow(result_table)
  # Number of missing values of estimated beta1 (e.g due to convergence)
  missing_gamma1 <- sum(is.na(result_table$est_E))
  
  # Bias
  gamma_1_hat <- mean(result_table[,"est_E"])
  IRR_E_hat <- mean(result_table[,"IRR_E"])
  bias_gamma_1 <- gamma_1_hat - true_gamma1 #absolute bias log scale
  se_gamma_1 <- sqrt(1/(n_sim-1)*sum((result_table[,"est_E"] - gamma_1_hat)^2)) #Empirical standard error
  bias_gamma_1_MCSE <- sqrt(1/n_sim)*se_gamma_1 #Monte Carlo standard error (MCSE) of absolute bias
  se_gamma_1_MCSE <- se_gamma_1/sqrt(2*(n_sim-1))
  
  percent_bias_g1 <- abs(bias_gamma_1)/true_gamma1
  bias_IRR_E <- IRR_E_hat - true_IRR_E #absolute bias IRR scale 
  
  gamma_2_hat <- mean(result_table[,"est_C"])
  bias_gamma_2 <- gamma_2_hat - true_gamma2
  percent_bias_g2 <- abs(bias_gamma_2)/true_gamma2
  
  # Coverage
  result_table$coverage_IRR_E <- with(result_table,
                                      IRR_E_CI_Lower <= true_IRR_E & IRR_E_CI_Upper >= true_IRR_E)
  coverage_irr_E <- mean(result_table$coverage_IRR_E)
  coverage_irr_E_MCSE <- sqrt(coverage_irr_E*(1-coverage_irr_E)/n_sim)
  
  result_table$coverage_IRR_C <- with(result_table,
                                      IRR_C_CI_Lower <= true_IRR_C & IRR_C_CI_Upper >= true_IRR_C)
  coverage_irr_C <- mean(result_table$coverage_IRR_C)  
  
  performance <- data.frame(missing_gamma1, gamma_1_hat, IRR_E_hat,
                            bias_gamma_1, bias_gamma_1_MCSE,
                            se_gamma_1, se_gamma_1_MCSE,
                            percent_bias_g1, bias_IRR_E, 
                            coverage_irr_E, coverage_irr_E_MCSE, 
                            bias_gamma_2, percent_bias_g2, coverage_irr_C)
  
  performance
}

