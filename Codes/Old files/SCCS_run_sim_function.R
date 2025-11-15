##########################################
### Project: Time-varying confounding  ###
### affected by past exposure in SCCS  ###
##########################################

###################################
### Codes for running           ###
### the full simulation         ###
###################################

# Output structure: within the Folder "Results", five subfolders are created for
# six analysis methods (ignore_C, exclude_C, censor_C, adjust_C, 
# adjust_C_obs_at_exp). Within each "Method" folder, there will be 288 .csv file 
# corresponding to 288 scenarios. Each scenario contains 1000 rows of results of 
# 1000 runs (est_E, se_E, IRR_E,...). The structure is displayed as follows:

#./Results/
#│
#├── ignore_C/
#│    ├── S1_pC0.0005_ConE0.25_IRRE1_winC10_UonY0.25.csv
#│    ├── ...
#│    └── scenario_288.csv
#│
#├── ...
#│
#└── adjust_C_obs_at_exp/
#│    ├── S1_pC0.0005_ConE0.25_IRRE1_winC10_UonY0.25.csv
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
# 2 - List of methods 
##############################

# Need to be among 'ignore_C', 'exclude_C', 'censor_C', 'adjust_C', 'adjust_C_obs_at_exp', and 'adjust_calendar'
methods <- c("ignore_C",              # ignore C in model
             "exclude_C",             # exclude subjects with C
             "censor_C",              # censor observation time at first C occurrence
             "adjust_C",              # adjust for C in SCCS model
             "adjust_C_obs_at_exp", # adjust for C and start observation at exposure
             "adjust_calendar")       # adjust for calendar month in SCCS model

##############################
# 3 - Seed for each run 
##############################

get_seeds <- function(n_sim, scenario_table){
  set.seed(2025, kind = "L'Ecuyer-CMRG", sample.kind = "Rejection")
  n_seed <- nrow(scenario_table)*n_sim # 288 scenarios and 1000 replications per scenario
  seed <- sample(1:1e9, 
                 size = n_seed, 
                 replace = FALSE)
}


############################## 
# 4 - Simulation study in layers
##############################

# --- Layer 1: Perform one run ------------------------------------------------

perform_one_run <- function(seed, rep, scen, methods, output_dir) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  # generate data
  tryCatch({
    # --- Generate data ---
    data <- generate_full_cohort(
      n = 100000, obs_time = 500, risk_window_E = 30,
      p_U = 1e-3, U_on_C = 5, U_on_Y = scen[['UonY']],     
      p_C = scen[['p_C']], C_on_C = 0.3, E_on_C = 0.2, risk_window_C = scen[['risk_win_C']],
      p_E = 3.2e-3, C_on_E = scen[['ConE']], E_on_E = 1e-6,
      baseline_Y = 2e-5, IRR_E = scen[['IRR_E']], IRR_C = 5
    )
  }, error = function(e) {
    log_error(e, stage = "Data generation", seed = seed, scen_name = scen[['scen_name']], rep = rep)
    return(NULL) 
  })
  
  
  # Ensure directories exist once
  invisible(lapply(methods, function(a) {
    create_directory(file.path(output_dir, a))
  }))
  
  # perform 6 analysis methods and output results 
  for (meth in methods) {

      file_path <- file.path(output_dir, meth, paste0(scen[['scen_name']], ".csv"))
      
      tryCatch({
      data_SCCS <- SCCS_reformat(data = data, risk_window_E = 30, 
                                 risk_window_C = scen[['risk_win_C']], method = meth)
      }, error = function(e) {
        log_error(e, stage = "Reformatting",
                  seed = seed, scen_name = scen[['scen_name']], rep = rep, method = meth)
      })
      
      tryCatch({      
      out <- analyse_sccs(data = data_SCCS, rep = rep, method = meth)
      }, error = function(e) {
        log_error(e, stage = "Analysis",
                  seed = seed, scen_name = scen[['scen_name']], rep = rep, method = meth)
      })
      
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

full_simulation <- function(scenario_table = scenarios, 
                            n_sim = n_sim, 
                            seeds = get_seeds(n_sim = n_sim, scenario_table = scenarios), 
                            methods = methods, 
                            output_dir = here("Results")) {
  
  create_directory(output_dir)
  for (s in seq_len(nrow(scenario_table))){
    invisible(
      sim_one_scenario(scen = scenario_table[s,], 
                       n_sim = n_sim, 
                       seeds = seeds, 
                       methods = methods, 
                       output_dir = output_dir)
    )
  }
}

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


#plan(multisession, workers = 4)
#n_sim = 4

#tic("test all scenarios")
#full_simulation(scenario_table = scenarios[30:nrow(scenarios),], 
#                n_sim = n_sim, 
#                seeds = get_seeds(n_sim, scenario_table = scenarios), 
#                methods = methods, 
#                output_dir = here("Results"))
#toc()

############################## 
# 5 - Summarise the simulation results
##############################

# Tables of all methods and scenarios

method_scen <- function(method_table = as.data.frame(methods), 
                        scenario_table = scenarios ){
  merged_table <- merge(method_table, scenario_table, by = NULL)
  return(merged_table)
}
  

# Function to read in the .csv simulation result files and summarise the bias

summarise_simulation_results <- function(method_scen = method_scen(), 
                                         results_dir = here("Results"),
                                         summary_dir = file.path(results_dir, "Summary"),
                                         summary_file_name = "Summary_all_scens") {
  create_directory(summary_dir)
  
  all_summaries <- list()
  
  # Read in .csv files
  for (i in seq_len(nrow(method_scen))) {
    method_i <- method_scen$methods[i]
    scen_file <- paste0(method_scen$scen_name[i], ".csv")
    
    file_path <- file.path(results_dir, method_i, scen_file)
    
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    sim_data <- read.csv(file_path)
    if (is.null(sim_data)) next
    
    # Quantify bias
    summary_row <- tryCatch(
      bias_quantification(true_IRR_E = method_scen$IRR_E[i],
                          true_IRR_C = 5,
                          result_table = sim_data),
      error = function(e) {
        warning(paste("Error applying bias_quantification to", file_path, ":", e$message))
        return(NULL)
      }
    )
    if (is.null(summary_row)) next
    
    # Combine with scenario information
    summary_row$methods <- method_i
    summary_row$scen_name <- method_scen$scen_name[i]

    
    
    all_summaries[[length(all_summaries) + 1]] <- summary_row
  }
  
  # Combine all summaries
  summary_table <- do.call(rbind, all_summaries)
  if (!"methods" %in% names(summary_table)) summary_table$methods <- NA
  if (!"scen_name" %in% names(summary_table)) summary_table$scen_name <- NA
  summary_table_merge <- merge(method_scen, summary_table, 
                               by = c("methods", "scen_name"),
                               all.x = TRUE)
  
  # Export file
  output_path <- file.path(summary_dir, paste0(summary_file_name,".csv"))
  write.csv(summary_table_merge, output_path, row.names = FALSE)
  
  return(summary_table_merge)
}

test <- summarise_simulation_results(method_scen = method_scen(scenario_table = scenarios[30:nrow(scenarios),]))
