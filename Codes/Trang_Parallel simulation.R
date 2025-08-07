#Load packages and data
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
pacman::p_load(
  foreach,     # foreach loop
  stats,
  tidyverse,
  SCCS,
  tictoc,      # Measure performance time
  progressr,   # Measure progress
  rio,         # Export file  
  here,
  data.table,
  magrittr,     # To use the pipe %>%
  doRNG,        # Reproducible parallel session
  future
) 

options(scipen = 999)

##########---Functions---####################################################

# Generate data --------------------------------

cohort_time_var_past_e_u <-function(n = 100000, obs_time = 500, risk_window = chosen_risk_window,
                                    p_U = 1e-3, U_on_C = 5, U_on_Y = 5,     
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
  
  U <- matrix(0L, nrow = n, ncol = obs_time)
  C <- matrix(0L, nrow = n, ncol = obs_time)
  E <- matrix(0L, nrow = n, ncol = obs_time)
  Y <- matrix(0L, nrow = n, ncol = obs_time)
  
  # Rolling sum of history of U, C and E
  
  roll_U_onC <- integer(n)
  roll_U_onY <- integer(n)
  roll_C_onC <- integer(n)
  roll_E_onC <- integer(n)
  roll_C_onE <- integer(n)
  roll_E_onE <- integer(n)
  roll_C_onY <- integer(n)
  roll_E_onY <- integer(n)
  
  
  for (t in 1:obs_time) {
    # Remove U, C or E events that fall outside the look-back window (28 days for C on E, C on Y, E on Y, and 91 days for U on C, U on Y, C on C and E on C)
    if (t > 28) {
      roll_C_onE <- roll_C_onE - C[, t- 28]
      roll_C_onY <- roll_C_onY - C[, t - 28]
      
    }
    if (t > risk_window) {
      roll_E_onY <- roll_E_onY - E[, t - risk_window]
    }
    if (t > 91) {
      roll_U_onC <- roll_U_onC - U[, t - 91]
      roll_U_onY <- roll_U_onY - U[, t - 91]                            
      roll_C_onC <- roll_C_onC - C[, t - 91]
      roll_E_onC <- roll_E_onC - E[, t - 91]
    }
    
    # Generate U
    prob_U <- plogis(log_p_U)
    U_t <- rbinom(n, 1, prob_U)
    U[, t] <- U_t
    
    # Update rolling of U
    roll_U_onC <- roll_U_onC + U_t
    roll_U_onY <- roll_U_onY + U_t
    
    # Generate C (depends on U[t-89, t], C[t-90, t-1] and E[t-90, t-1])
    recent_U_onC <- roll_U_onC > 0 # Check recent history of U  
    recent_C_onC <- roll_C_onC > 0 # Check recent history of C
    recent_E_onC <- roll_E_onC > 0 # Check recent history of E
    prob_C <- plogis(log_p_C + log_U_on_C*recent_U_onC + log_C_on_C*recent_C_onC + log_E_on_C*recent_E_onC)
    C_t <- rbinom(n, 1, prob_C)
    C[, t] <- C_t
    
    # Update rolling of C 
    roll_C_onC <- roll_C_onC + C_t
    roll_C_onE <- roll_C_onE + C_t
    roll_C_onY <- roll_C_onY + C_t
    
    # Generate E (depends on C[t-27, t] and E[1, t-1])
    recent_C_onE <- roll_C_onE > 0
    recent_E_onE <- roll_E_onE > 0
    prob_E <- plogis(log_p_E + log_C_on_E*recent_C_onE + log_E_on_E*recent_E_onE)
    E_t <- rbinom(n, 1, prob_E)
    E[, t] <- E_t
    
    # Update rolling of E
    roll_E_onE <- roll_E_onE + E_t
    roll_E_onC <- roll_E_onC + E_t
    roll_E_onY <- roll_E_onY + E_t
    
    
    # Generate Y (depends on U[t-89, t] C[t-27, t], E[t-risk_window+1, t])
    recent_U_onY <- roll_U_onY > 0
    recent_C_onY <- roll_C_onY > 0
    recent_E_onY <- roll_E_onY > 0
    
    prob_Y <- plogis(log_baseline_Y + log_U_on_Y*recent_U_onY + log_IRR_E*recent_E_onY + log_IRR_C*recent_C_onY)
    Y[, t] <- rbinom(n, 1, prob_Y)
    
  }
  
  # Reshape to long format
  data_long <- data.table(
    id = rep(1:n, each = obs_time),
    day = rep(1:obs_time, times = n),
    U = as.vector(t(U)), #transpose vector
    C = as.vector(t(C)),
    E = as.vector(t(E)),
    Y = as.vector(t(Y))
  )
  
  return(data_long)
}

# Reshape data --------------------------------------------------------------

SCCS_reformat_confound <- function(data, risk_window = chosen_risk_window){
  Y_data <- data %>%
    .[Y ==1 , .(id, Y_day = day)] %>%
    .[, .row := seq_len(.N), by = id]
  E_data <- data %>%
    .[E ==1 , .(id, E_day = day)] %>%
    .[, .row := seq_len(.N), by = id]
  C_data <- data %>%
    .[C ==1 , .(id, C_day = day)] %>%
    .[, .row := seq_len(.N), by = id]
  
  merge_E_C <- merge(E_data, C_data, by = c("id", ".row"), all = TRUE)
  merge_Y_E_C <- merge(Y_data, merge_E_C, by = c("id"), all = TRUE)
  
  merge_Y_E_C <- merge_Y_E_C %>%
    .[!is.na(Y_day)] %>%
    .[ , `:=`(obs_start = 1, 
              obs_end = 500,
              E_end = E_day + risk_window - 1,
              C_end = C_day + 28 - 1)]
  
  merge_Y_E_C <- as.data.frame(merge_Y_E_C)
  
  return(merge_Y_E_C)
}

# Analyse data --------------------------------------------------------------

analyse_sccs_confound <- function(data, rep)
{
  # Fit SCCS model
  model <- standardsccs(event ~ E_day + C_day, 
                        indiv= id,         # subject ID
                        astart = obs_start,# start of observation period 
                        aend = obs_end,    # end of observation period 
                        aevent = Y_day,    # event time
                        adrug = cbind(E_day, C_day), # start of exposure
                        aedrug = cbind(E_end, C_end),  # end of risk period
                        expogrp = list(0, 0),       # start of risk period counted from 'adrug'
                        data = data)
  
  # Extract results from SCCS model
  est_E <- coef(model)[1,1]
  se_E <- coef(model)[1,3]
  IRR_E <- exp(est_E)
  IRR_E_CI_Lower <- model$conf.int[1,3]
  IRR_E_CI_Upper <- model$conf.int[1,4]
  est_C <- coef(model)[2,1]
  IRR_C <- coef(model)[2,2]
  IRR_C_CI_Lower <- model$conf.int[2,3]
  IRR_C_CI_Upper <- model$conf.int[2,4]
  n_event <- model$nevent
  res <- data.frame(rep = rep, 
                    est_E, se_E, IRR_E, IRR_E_CI_Lower, IRR_E_CI_Upper, 
                    est_C, IRR_C, IRR_C_CI_Lower, IRR_C_CI_Upper, 
                    n_event, 
                    row.names = NULL)
  return(res)
}

##################--Run the simulation in parallel--##############################

set.seed(1997) 
registerDoRNG(1997) 

# Set up parallel sessions
plan(multisession, workers = 10)

n_sim <- 1000
chosen_risk_window <- 28

tic("Parallel simulation")
results_list <- foreach(i = 1:n_sim, 
                        .options.future = list(packages = c("extraDistr", "dplyr", "SCCS", "data.table"),
                                               seed = TRUE),
                        .combine = rbind) %dofuture% {
                          data <- cohort_time_var_past_e_u(n = 100000, obs_time = 500, risk_window = chosen_risk_window,
                                                           p_U = 1e-3, U_on_C = 5, U_on_Y = 5,     
                                                           p_C = 1e-3, C_on_C = 0.3, E_on_C = 0.2,
                                                           p_E = 3.2e-3, C_on_E = 0.25, E_on_E = 0.001,
                                                           baseline_Y = 2e-5, IRR_E = 2, IRR_C = 5)
                          data_SCCS <- SCCS_reformat_confound(data)
                          result<-analyse_sccs_confound(data= data_SCCS, rep = i)
                          result
                        }

toc() # measure run time 