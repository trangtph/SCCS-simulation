##########################################
### Project: Time-varying confounding  ###
### affected by past exposure in SCCS  ###
##########################################

###################################
### Functions for               ###
### simulating and analyzing    ###
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

# ------------------------------------------------------------------------------
# List of methods --------------------------------------------------------------
# ------------------------------------------------------------------------------
#methods <- c("ignore_C",              # ignore C in model
#             "exclude_C",             # exclude subjects with C
#             "censor_C",              # censor observation time at first C occurrence
#             "adjust_C",              # adjust for C in SCCS model
#             "adjust_C_obs_at_exp", # adjust for C and start observation at exposure
#             "adjust_calendar")       # adjust for calendar month in SCCS model

# ------------------------------------------------------------------------------
# Function: generate cohort data -----------------------------------------------
# ------------------------------------------------------------------------------
# Four variables: exposure E, outcome Y, confounder C of E and Y that is affected
# by past E, confounder U of C and Y

generate_full_cohort <-function(n = 100000, obs_time = 500, risk_window_E = 30,
                                    p_U = 1e-3, U_on_C = 5, U_on_Y = 15,     
                                    p_C = 1e-3, C_on_C = 0.3, E_on_C = 0.2, risk_window_C = chosen_risk_window_C,
                                    p_E = 3.2e-3, C_on_E = 0.25, E_on_E = 1e-6,
                                    baseline_Y = 2e-5, IRR_E = 2, IRR_C = 5)
{
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
    # (30 days for C on E, /risk window_C/ for C on Y, /risk window_E/ for E on Y,
    # 90 days for U on C, U on Y, and 91 days for C on C and E on C)
    if (t > 30) {
      roll_C_onE <- roll_C_onE - C[, t- 30]
    }
    if (t > risk_window_C) {
      roll_C_onY <- roll_C_onY - C[, t - risk_window_C]
    }
    if (t > risk_window_E) {
      roll_E_onY <- roll_E_onY - E[, t - risk_window_E]
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
    U_t <- stats::rbinom(n, 1, p_U)
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
    
    
    # Generate Y (depends on U[t-89, t] C[t-risk_windowC+1, t], E[t-risk_windowE+1, t])
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
# Only keep data on date of occurrence of exposure, outcome and covariates
# Method ignore_C, exclude_C, censor_C, adjust_C, C_in_model_obs_at_exp, adjust_calendar

SCCS_reformat <- function(data, risk_window_E = 30, risk_window_C = chosen_risk_window_C, method){
  Y_data <- data[Y == 1, .(id, Y_day = day)] #Filter rows with Y==1, select id and day columns
  Y_data[, .row := seq_len(.N), by = id] #For each individual, create event index
  
  E_data <- data[E == 1, .(id, E_day = day)]
  E_data[, .row := seq_len(.N), by = id]
  E_data[, E_end := pmin(500, E_day + risk_window_E - 1)]
  E_data <- E_data[, .SD[1], by = .(id, E_end)]
  
  C_data <- data[C == 1, .(id, C_day = day)]
  C_data[, .row := seq_len(.N), by = id]
  C_data[, C_end := pmin(500, C_day + risk_window_C - 1)]
  C_data <- C_data[, .SD[1], by = .(id, C_end)] #if there are >1 records within 1 id with the same C_end, keep the smallest C_day
  # rationale for the above step is because if there are 2 risk windows from the same exposure in which one is totally inside the other,
  # it will cause error when fitting standardSCCS()
  
  if(method %in% c("ignore_C", "adjust_calendar")){
    merge_Y_E <- data.table::merge.data.table(Y_data, E_data, by = "id", all = TRUE)
    merge_Y_E <- merge_Y_E[!is.na(Y_day)] # Keep only rows with observed events
    
    # Add observation window and end of risk periods
    merge_Y_E[, `:=`(
      obs_start = 1,
      obs_end = 500
    )]
    return(as.data.frame(merge_Y_E))
    
  }
  
  else if(method == "exclude_C"){ #Exclude subjects with C
    subject_C <- unique(C_data$id) #Create a list of subjects with C
    
    merge_Y_E <- data.table::merge.data.table(Y_data, E_data, by = "id", all = TRUE)
    merge_Y_E <- merge_Y_E[!is.na(Y_day)]
    
    merge_Y_E[, `:=`(
      obs_start = 1,
      obs_end = 500
    )]    
    
    merge_Y_E_noC <- merge_Y_E[!id %in% subject_C] #Exclude subjects with C
    return(as.data.frame(merge_Y_E_noC))
  }
  
  else if(method == "censor_C"){ # Censor subject on the day of first C
    C_data_1st <- C_data[, .SD[1], by = id] # Retain the first C of each subject
    C_data_1st <- C_data_1st[,.(id, C_day_censor = C_day)]
    
    merge_Y_E <- data.table::merge.data.table(Y_data, E_data, by = "id", all = TRUE)
    merge_Y_E <- merge_Y_E[!is.na(Y_day)]
    
    merge_Y_E[, `:=`(
      obs_start = 1,
      obs_end = 500
    )]
    
    merge_Y_E_Ccensored <- merge_Y_E[C_data_1st, on = .(id), obs_end := C_day_censor-1]
    merge_Y_E_Ccensored[, E_day := fifelse(E_day > obs_end, NA, E_day)]
    merge_Y_E_Ccensored[, E_end := fifelse(is.na(E_day), NA, E_end)]

    return(as.data.frame(merge_Y_E_Ccensored))
  }
  
  else if(method == "adjust_C"){ # Adjust for C in the model
    
    merge_E_C <- data.table::merge.data.table(E_data, C_data, by = c("id", ".row"), all = TRUE)
    merge_Y_E_C <- data.table::merge.data.table(Y_data, merge_E_C, by = "id", all = TRUE)
    merge_Y_E_C <- merge_Y_E_C[!is.na(Y_day)]
    
    merge_Y_E_C[, `:=`(
      obs_start = 1,
      obs_end = 500
    )]
    
    return(as.data.frame(merge_Y_E_C))
  }
  
  else if(method == "adjust_C_obs_at_exp"){ # Start observation at exposure
    
    E_data_1st <- E_data[, .SD[1], by = id] # Retain the first E of each subject
    E_data_1st <- E_data_1st[,.(id, E_day_censor = E_day)]
    
    merge_E_C <- data.table::merge.data.table(E_data, C_data, by = c("id", ".row"), all = TRUE)
    merge_Y_E_C <- data.table::merge.data.table(Y_data, merge_E_C, by = "id", all = TRUE)
    merge_Y_E_C <- merge_Y_E_C[!is.na(Y_day)]
    
    merge_Y_E_C[, `:=`(
      obs_start = 1,
      obs_end = 500
    )]

    merge_Y_E_C_left_censored <- merge_Y_E_C[E_data_1st, on = .(id), obs_start := E_day_censor] #Start observation on the day of 1st E
    merge_Y_E_C_left_censored[, C_day := case_when(!is.na(C_day) & C_day<obs_start-risk_window_C +1 ~ NA , #Left-censor C on obs_start-risk_window_C +1
                                                   TRUE ~ C_day)]
    merge_Y_E_C_left_censored[, `:=`(
      E_end = fifelse(is.na(E_day), NA, E_end),
      C_end = fifelse(is.na(C_day), NA, C_end)
    )]
    
    merge_Y_E_C_left_censored[, C_day := fifelse(!is.na(C_day) & C_day <= obs_start, obs_start, C_day)] #Left-censor start of C on obs_start
    
    
    return(as.data.frame(merge_Y_E_C_left_censored))
  }
  
  else {
    stop("Error: Method should be among ignore_C, exclude_C, censor_C, adjust_C, adjust_C_obs_at_exp, and adjust_calendar")
  }
    
  }

# Function: analyse data ------------------------------------------------------
# Three models: ignore C (method ignore_C, exclude_C, censor_C), 
#               adjust for C (method adjust_C, adjust_C_obs_at_exp), 
#               and adjust for calendar time (method adjust_calendar)

analyse_sccs <- function(data, rep, method)
{
  if(method %in% c("adjust_C", "adjust_C_obs_at_exp")){
    # Adjust for C
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
  
  
  else if(method == "adjust_calendar"){
    # Adjust for calendar time
    model <- SCCS::standardsccs(event ~ E_day + age, 
                                indiv= id,      # subject ID
                                astart = obs_start,# start of observation period 
                                aend = obs_end,    # end of observation period 
                                aevent = Y_day, # event time
                                adrug = E_day, # start of exposure
                                aedrug = E_end,  # end of risk period
                                expogrp = 0,       # start of risk period counted from 'adrug'
                                agegrp = seq(31, 500, by = 30),
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
  else{
    # Model with only effect of E
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
  

  # Non-convergence issues
  non_converge_gamma_1 <- sum(is.na(result_table$est_E), nrow(result_table[result_table$est_E <= -10 | result_table$est_E >= 10,]), na.rm = TRUE) 
  result_table <- result_table[result_table$est_E > -10 & result_table$est_E < 10,]
  n_sim <- nrow(result_table)
  
  #mean number of events
  mean_nr_event <- mean(result_table$n_event, na.rm = TRUE)
  # Estimates
  gamma_1_hat <- mean(result_table[,"est_E"], na.rm = TRUE)
  IRR_E_hat <- mean(result_table[,"IRR_E"], na.rm = TRUE)
  
  # Variance 
  se_gamma_1 <- sqrt(1/(n_sim-1)*sum((result_table[,"est_E"] - gamma_1_hat)^2, na.rm = TRUE)) #Empirical standard error
  se_gamma_1_MCSE <- se_gamma_1/sqrt(2*(n_sim-1)) # Monte Carlo standard error (MCSE) of empirical SE
  se_gamma_1_low_CI <- se_gamma_1 - 1.96*se_gamma_1_MCSE
  se_gamma_1_up_CI <- se_gamma_1 + 1.96*se_gamma_1_MCSE
  
  mod_se_gamma_1 <- sqrt(mean((result_table[,"se_E"])^2, na.rm = TRUE)) # Model-based SE, target empirical SE
  
  # Bias
  bias_gamma_1 <- mean(result_table[,"est_E"]-true_gamma1, na.rm = TRUE) #absolute bias log scale
  bias_gamma_1_MCSE <- sqrt(1/n_sim)*se_gamma_1 #  MCSE of absolute bias
  bias_gamma_1_low_CI <- bias_gamma_1 - 1.96*bias_gamma_1_MCSE
  bias_gamma_1_up_CI <- bias_gamma_1 + 1.96*bias_gamma_1_MCSE
  
  gamma_1_hat_low_CI <-  gamma_1_hat - 1.96*bias_gamma_1_MCSE
  gamma_1_hat_up_CI <-  gamma_1_hat + 1.96*bias_gamma_1_MCSE
  
  relative_bias_g1 <- mean((result_table[,"est_E"]-true_gamma1)/true_gamma1, na.rm = TRUE) # relative bias log scale
  relative_bias_g1_MCSE <- sqrt(1/(n_sim*(n_sim-1))*sum(((result_table[,"est_E"]-true_gamma1)/true_gamma1-relative_bias_g1)^2, na.rm = TRUE))
  relative_bias_g1_low_CI <- relative_bias_g1 - 1.96*relative_bias_g1_MCSE
  relative_bias_g1_up_CI <- relative_bias_g1 + 1.96*relative_bias_g1_MCSE
  
  bias_IRR_E <- mean(result_table[,"IRR_E"]-true_IRR_E, na.rm = TRUE) #absolute bias IRR scale 
  bias_IRR_E_MCSE <- sqrt(1/((n_sim-1)*n_sim)*sum((result_table[,"IRR_E"] - IRR_E_hat)^2, na.rm = TRUE))
  bias_IRR_E_low_CI <- bias_IRR_E - 1.96*bias_IRR_E_MCSE
  bias_IRR_E_up_CI <- bias_IRR_E + 1.96*bias_IRR_E_MCSE
  IRR_E_hat_low_CI <- IRR_E_hat - 1.96*bias_IRR_E_MCSE
  IRR_E_hat_up_CI <- IRR_E_hat + 1.96*bias_IRR_E_MCSE
  
  gamma_2_hat <- mean(result_table[,"est_C"], na.rm = TRUE)
  bias_gamma_2 <- gamma_2_hat - true_gamma2
  percent_bias_g2 <- abs(bias_gamma_2)/true_gamma2
  
  #MSE
  mse_g1 <- mean((result_table[,"est_E"]-true_gamma1)^2, na.rm = TRUE)
  mse_MCSE <- sqrt(1/(n_sim-1)*mean(((result_table[,"est_E"]-true_gamma1)^2 -mse_g1)^2, na.rm = TRUE))
  mse_g1_low_CI <- mse_g1 - 1.96*mse_MCSE
  mse_g1_up_CI <- mse_g1 + 1.96*mse_MCSE
  
  # Coverage
  result_table$coverage_IRR_E <- with(result_table,
                                      IRR_E_CI_Lower <= true_IRR_E & IRR_E_CI_Upper >= true_IRR_E)
  coverage_irr_E <- mean(result_table$coverage_IRR_E)
  coverage_irr_E_MCSE <- sqrt(coverage_irr_E*(1-coverage_irr_E)/n_sim)
  coverage_irr_E_low_CI <- coverage_irr_E - 1.96*coverage_irr_E_MCSE
  coverage_irr_E_up_CI <- coverage_irr_E + 1.96*coverage_irr_E_MCSE
  
  result_table$coverage_IRR_C <- with(result_table,
                                      IRR_C_CI_Lower <= true_IRR_C & IRR_C_CI_Upper >= true_IRR_C)
  coverage_irr_C <- mean(result_table$coverage_IRR_C)  
  
  # Summary table
  performance <- data.frame(mean_nr_event, non_converge_gamma_1, 
                            gamma_1_hat, gamma_1_hat_low_CI, gamma_1_hat_up_CI,
                            IRR_E_hat, IRR_E_hat_low_CI, IRR_E_hat_up_CI,
                            bias_gamma_1, bias_gamma_1_MCSE, bias_gamma_1_low_CI, bias_gamma_1_up_CI,
                            mod_se_gamma_1, se_gamma_1, se_gamma_1_MCSE, se_gamma_1_low_CI, se_gamma_1_up_CI,
                            relative_bias_g1, relative_bias_g1_MCSE, relative_bias_g1_low_CI, relative_bias_g1_up_CI,
                            bias_IRR_E, bias_IRR_E_MCSE, bias_IRR_E_low_CI, bias_IRR_E_up_CI,
                            mse_g1, mse_MCSE, mse_g1_low_CI, mse_g1_up_CI,
                            coverage_irr_E, coverage_irr_E_MCSE, coverage_irr_E_low_CI, coverage_irr_E_up_CI,
                            bias_gamma_2, percent_bias_g2, coverage_irr_C)
  
  return(performance)
}