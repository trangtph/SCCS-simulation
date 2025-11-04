###################################
### Codes for running           ###
### the full simulation         ###
###################################

# Output structure: within the Folder "Results", six subfolders are created for
# six analysis methods (ignore_C, exclude_C, censor_C, adjust_C, 
# C_in_model_obs_at_exp, adjust_calendar). Within each "Method" folder, there will be 288 .csv file 
# corresponding to 288 scenarios. Each scenario contains 1000 rows of results of 
# 1000 runs (est_E, se_E, IRR_E,...). The structure is displayed as follows:

#./Results/
#│
#├── method1_ignore_C/
#│    ├── scenario_1.csv
#│    ├── ...
#│    └── scenario_288.csv
#│
#├── ...
#│
#└── method6_calendar/
#│    ├── scenario_1.csv
#│    ├── ...
#│    └── scenario_288.csv

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

source(here("Codes", "SCCS_sim_functions.R"))
source(here("Codes", "SCCS_helper_functions.R"))

options(scipen = 999)

##############################
# 1 - Table of 288 scenarios 
##############################

# Parameter grids
p_C_vals    <- c(5e-4, 1e-3, 5e-3)
C_on_E_vals <- c(0.25, 1, 2, 5)
IRR_E_vals  <- c(1, 2)
chosen_risk_window_C <- c(10, 30, 60)
U_on_Y_vals <- c(0.25, 1, 2, 15)

scenarios <- expand.grid(p_C = p_C_vals,
                         ConE = C_on_E_vals,
                         IRR_E = IRR_E_vals,
                         risk_win_C = chosen_risk_window_C,
                         UonY = U_on_Y_vals)

scenarios$scen_id <- seq(1:nrow(scenarios))
scenarios$scen_name <- with(scenarios, 
                                paste0("S",scen_id,"_pC", p_C, "_ConE", ConE,
                                       "_IRRE", IRR_E, "_winC",risk_win_C, "_UonY", UonY))

##############################
# 2 - Seed for each run 
##############################

get_seeds <- function(n_sim){
  set.seed(2025)
  n_seed <- nrow(scenarios)*n_sim # 288 scenarios and 1000 replications per scenario
  seed <- sample(1:1e9, 
                 size = n_seed, 
                 replace = FALSE)
}

############################## 
# 3 - Simulation study in layers
##############################

# --- Layer 1: Perform one run ------------------------------------------------

perform_one_run <- function(seed, rep, scen, methods, output_dir) {
  set.seed(seed)
  # generate data
  data <- generate_full_cohort(n = 100000, obs_time = 500, risk_window_E = 30,
                               p_U = 1e-3, U_on_C = 5, U_on_Y = scen[['UonY']],     
                               p_C = scen[['p_C']], C_on_C = 0.3, E_on_C = 0.2, risk_window_C = scen[['risk_win_C']],
                               p_E = 3.2e-3, C_on_E = scen[['ConE']], E_on_E = 1e-6,
                               baseline_Y = 2e-5, IRR_E = scen[['IRR_E']], IRR_C = 5)
  
  # Ensure directories exist once
  invisible(lapply(methods, function(a) {
    create_directory(file.path(output_dir, a))
  }))
  
  # perform 6 analysis methods and output results 
  for (meth in methods) {
    file_path <- file.path(output_dir, meth, paste0(scen[['scen_name']], ".csv"))
    
    data_SCCS <- SCCS_reformat(data = data, risk_window_E = 30, 
                               risk_window_C = scen[['risk_win_C']], method = meth)
    out <- analyse_sccs(data = data_SCCS, rep = rep, method = meth)
    out <- as.data.frame(cbind(out, seed = seed))
    append_to_csv(out, file_path)
  }
}

# --- Layer 2: Repeat 'Perform one run' n_sim times, for each scenario ---------
# (Parallelized)

sim_one_scenario <- function(scen, n_sim, seeds, methods, output_dir){
  
  message(paste("Running scenario:", scen[['scen_id']], scen[['scen_name']]))
  
  scen_num <- as.numeric(scen['scen_id'])

  foreach(i = 1:n_sim,
          .options.future = list(packages = c("extraDistr","dplyr","SCCS","data.table","RhpcBLASctl"),
                                            seed = TRUE)) %dofuture% {
                                              
           # To avoid issues with native libraries, try keeping native libraries thread safe
           data.table::setDTthreads(1L)
           RhpcBLASctl::blas_set_num_threads(1L)
           RhpcBLASctl::omp_set_num_threads(1L)
           
           seed = seeds[(n_sim *(scen_num -1) + i)] #S1 gets seed nr 1-1000, S2 get seed nr 1001-2000 and so on
           perform_one_run(seed = seed, rep = i, scen = scen, methods = methods, output_dir = output_dir)
           }
  }

# --- Layer 3: Perform the full simulation ------------------------------------
# For each scenario, generate n_sim datasets and perform 6 analysis methods

full_simulation <- function(scenario_list = scenarios, 
                            n_sim = n_sim, 
                            seeds = get_seeds(n_sim), 
                            methods = methods, 
                            output_dir = here("Results")) {
  
  create_directory(output_dir)
  for (s in seq_len(nrow(scenario_list))){
    invisible(
      sim_one_scenario(scen = scenario_list[s,], n_sim = n_sim, seeds = seeds, methods = methods, output_dir = output_dir)
    )
  }
}

############################## 
# 4 - Test the simulation
##############################

perform_one_run(seed = 1357, rep = 1, scen = scenarios[1,], methods = methods, output_dir = here("Results"))

plan(multisession, workers = 2)
sim_one_scenario(scen = scenarios[1,], n_sim=3, seeds = get_seeds(n_sim=3), methods = methods, output_dir = here("Results"))

scen_test <- scenarios[1:4,]

plan(multisession, workers = 4)
n_sim = 4
full_simulation(scenario_list = scenarios, 
                n_sim = n_sim, 
                seeds = get_seeds(n_sim), 
                methods = methods, 
                output_dir = here("Results"))


