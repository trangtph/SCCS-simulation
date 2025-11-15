##########################################
### Project: Time-varying confounding  ###
### affected by past exposure in SCCS  ###
##########################################

###################################
### Codes for running           ###
### the full simulation         ###
###################################

##############################
# 0 - Load libraries and source code 
##############################

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
  progressr,   # Measure progress
  rio,         # Export file  
  here,
  data.table,
  magrittr,     # To use the pipe %>%
  doRNG,        # Reproducible parallel session
  RhpcBLASctl   # Control threads in parallel session
) 

source(here("Codes", "SCCS_main_functions.R"))
source(here("Codes", "SCCS_sim_workflow_functions.R"))
source(here("Codes", "SCCS_helper_functions.R"))

options(scipen = 999)

##############################
# 1 - Table of all scenarios 
##############################

scenarios <- scenario_table(p_C_vals = c(5e-4, 1e-3, 5e-3), 
                             C_on_E_vals = c(0.25, 1, 2, 5),
                             IRR_E_vals = c(1, 2), 
                             chosen_risk_window_C = c(10, 30, 60),
                             U_on_Y_vals = c(0.25, 1, 2, 15))

##############################
# 2 - List of methods 
##############################

# Need to be among 'ignore_C', 'exclude_C', 'censor_C', 'adjust_C', 'adjust_C_obs_at_exp', and 'adjust_calendar'
methods <- c("ignore_C",              # ignore C in model
             "exclude_C",             # exclude subjects with C
             "censor_C",              # censor observation time at first C occurrence
             "adjust_C",              # adjust for C in SCCS model
             "adjust_C_obs_at_exp") # adjust for C and start observation at exposure

##############################
# 3 - Generate data and apply analysis methods 
##############################
plan(multisession, workers = 4)
n_sim <- 1000
full_simulation(scenario_table = scenarios, 
                n_sim = n_sim, 
                seeds = get_seeds(n_sim, scenario_table = scenarios), 
                methods = methods, 
                output_dir = here("Results"))


##############################
# 4 - Quantify bias and summarise results
##############################
bias_results <- summarise_simulation_results(method_scen = method_scen(scenario_table = scenarios),
                                             results_dir = here("Results"),
                                             summary_dir = file.path(results_dir, "Summary"),
                                             summary_file_name = "Summary_all_scens")



############################## 
# Test the simulation
##############################

#perform_one_run(seed = 1357, rep = 1, scen = scenarios[1,], methods = methods, output_dir = here("Results"))

#plan(multisession, workers = 4)
#sim_one_scenario(scen = scenarios[1,], 
#                 n_sim=3, 
#                 seeds = get_seeds(n_sim=3, scenario_table = scenarios), 
#                 methods = methods,
#                 output_dir = here("Results"))




scenarios <- scenario_table(p_C_vals = c(5e-3), 
                            C_on_E_vals = c(5),
                            IRR_E_vals = c(2), 
                            chosen_risk_window_C = c(60),
                            U_on_Y_vals = c(2, 15))

plan(multisession, workers = 4)
n_sim <- 5
full_simulation(scenario_table = scenarios, 
                n_sim = n_sim, 
                seeds = get_seeds(n_sim, scenario_table = scenarios), 
                methods = methods, 
                output_dir = here("Results"))
bias_results <- summarise_simulation_results(method_scen = method_scen(scenario_table = scenarios))

