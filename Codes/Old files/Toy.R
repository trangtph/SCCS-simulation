scen <- scenarios[3,]
scen_num <- as.numeric(scen['scen_id'])
seed_val <- get_seeds(n_sim=4)
seed <- seed_val[(4 *(scen_num -1) + 2)] # run #2 (no error)
set.seed(seed)

data <- generate_full_cohort(n = 100000, obs_time = 500, risk_window_E = 30,
                              p_U = 1e-3, U_on_C = 5, U_on_Y = scen[['UonY']],     
                              p_C = scen[['p_C']], C_on_C = 0.3, E_on_C = 0.2, risk_window_C = scen[['risk_win_C']],
                              p_E = 3.2e-3, C_on_E = scen[['ConE']], E_on_E = 1e-6,
                              baseline_Y = 2e-5, IRR_E = scen[['IRR_E']], IRR_C = 5)
data_SCCS <- SCCS_reformat(data = data, risk_window_E = 30, 
                            risk_window_C = scen[['risk_win_C']], method = "adjust_C_obs_at_exp")
analyse_sccs(data = data_SCCS, rep = 1, method = "adjust_C_obs_at_exp")
# ------
seed <- seed_val[(4 *(scen_num -1) + 1)] # run #1 (subscript out of bounds error)
set.seed(seed)

data2 <- generate_full_cohort(n = 100000, obs_time = 500, risk_window_E = 30,
                             p_U = 1e-3, U_on_C = 5, U_on_Y = scen[['UonY']],     
                             p_C = scen[['p_C']], C_on_C = 0.3, E_on_C = 0.2, risk_window_C = scen[['risk_win_C']],
                             p_E = 3.2e-3, C_on_E = scen[['ConE']], E_on_E = 1e-6,
                             baseline_Y = 2e-5, IRR_E = scen[['IRR_E']], IRR_C = 5)

data_SCCS2 <- SCCS_reformat(data = data2, risk_window_E = 30, 
                           risk_window_C = scen[['risk_win_C']], method = "adjust_C_obs_at_exp")
analyse_sccs(data = data_SCCS2, rep = 1, method = "adjust_C_obs_at_exp")

# Reshape data function
data <- data2
Y_data <- data[Y == 1, .(id, Y_day = day)] #Filter rows with Y==1, select id and day columns
Y_data[, .row := seq_len(.N), by = id] #For each individual, create event index

E_data <- data[E == 1, .(id, E_day = day)]
E_data[, .row := seq_len(.N), by = id]

C_data <- data[C == 1, .(id, C_day = day)]
C_data[, .row := seq_len(.N), by = id]

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

merge_Y_E_C_left_censored[, `:=`(
  E_end = fifelse(is.na(E_day), NA, pmin(E_day + risk_window_E - 1, obs_end)),
  C_end = fifelse(is.na(C_day), NA, pmin(C_day + risk_window_C - 1, obs_end))
)]

data_SCCS4 <- as.data.frame(merge_Y_E_C_left_censored)

analyse_sccs(data = data_SCCS3, rep = 1, method = "adjust_C_obs_at_exp")

SCCS::standardsccs(event ~ E_day + C_day, 
                   indiv= id,         # subject ID
                   astart = obs_start,# start of observation period 
                   aend = obs_end,    # end of observation period 
                   aevent = Y_day,    # event time
                   adrug = cbind(E_day, C_day), # start of exposure
                   aedrug = cbind(E_end, C_end),  # end of risk period
                   expogrp = list(0, 0),       # start of risk period counted from 'adrug'
                   data = data_SCCS4)

bad <- data_SCCS %>% filter(
  (E_end < E_day) |
    (C_end < C_day) |
    (E_day < obs_start) |
    (E_end > obs_end) |
    (C_day < obs_start) |
    (C_end > obs_end) |
    (obs_end <= obs_start))
nrow(bad)

summary(data_SCCS)
summary(data_SCCS3)

any(data_SCCS$E_end <= data_SCCS$E_day)
any(data_SCCS$C_end <= data_SCCS$C_day)

any(data_SCCS3$Y_day > data_SCCS3$obs_end)

# -------------------------------
scen <- scenarios[30,]
scen_num <- as.numeric(scen['scen_id'])
seed_val <- get_seeds(n_sim=4)
seed <- seed_val[(4 *(scen_num -1) + 4)] # run #1 
set.seed(seed)

data <- generate_full_cohort(n = 100000, obs_time = 500, risk_window_E = 30,
                             p_U = 1e-3, U_on_C = 5, U_on_Y = scen[['UonY']],     
                             p_C = scen[['p_C']], C_on_C = 0.3, E_on_C = 0.2, risk_window_C = scen[['risk_win_C']],
                             p_E = 3.2e-3, C_on_E = scen[['ConE']], E_on_E = 1e-6,
                             baseline_Y = 2e-5, IRR_E = scen[['IRR_E']], IRR_C = 5)
data_SCCS <- SCCS_reformat(data = data, risk_window_E = 30, 
                           risk_window_C = scen[['risk_win_C']], method = "adjust_C_obs_at_exp")
analyse_sccs(data = data_SCCS, rep = 1, method = "adjust_C_obs_at_exp")