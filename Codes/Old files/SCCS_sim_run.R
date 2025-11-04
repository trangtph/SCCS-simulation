##########################################################
### Simulating and analysing data.                     ###
##########################################################

#Load packages and functions -------------------------------------------------------
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

source(here("SCCS_Simulation", "SCCS_sim_functions.R"))

options(scipen = 999)


# Approach 1: Ignoring C in SCCS model -----------------------------------------
#################################################################################
n_sim <- 1000   

# Parameter grids
p_C_vals    <- c(5e-4, 1e-3, 5e-3)
C_on_E_vals <- c(0.25, 1, 5)
IRR_C_vals  <- c(2, 5, 15)
IRR_E_vals  <- c(1, 2, 5)
chosen_risk_window <- c(10, 50)

scenarios <- expand.grid(p_C = p_C_vals,
                         IRR_C = IRR_C_vals,
                         IRR_E = IRR_E_vals,
                         ConE = C_on_E_vals,
                         risk_win = chosen_risk_window)

scenarios$scenario_id <- seq(1:nrow(scenarios))
scenarios$scenario_name <- with(scenarios, 
                                paste0("Scen_",scenario_id,"_pC", p_C, "_ConE", ConE,
                                       "_IRRC", IRR_C, "_IRRE",
                                       IRR_E, "_",risk_win, "d"))

#data.table::fwrite(scenarios, here("SCCS_Simulation", "Scenarios_ignore_C.csv"))

# Set the directory to export files
out_directory <- here::here("SCCS_Simulation", "Results", "Ignore_C")

## Run the first 15 scenarios ---------------------------------------------------

plan(multisession, workers = 40)

registerDoRNG()
set.seed(1010)

tic("Sim: bias when ignoring C")

scen_start <- 1
scen_stop <- 15
results_all <- foreach(s = scen_start:scen_stop,
                       .combine = rbind) %do% {
      p_C   <- scenarios$p_C[s]
      C_on_E <- scenarios$ConE[s]
      IRR_C <- scenarios$IRR_C[s]
      IRR_E <- scenarios$IRR_E[s]
      risk_win <- scenarios$risk_win[s]
      scen_name <- scenarios$scenario_name[s]
                         
      # Run simulation for one scenario
      sim_res <- foreach(i = 1:n_sim,
                        .combine = rbind,
                        .options.future = list(packages = c("extraDistr","dplyr","SCCS","data.table","RhpcBLASctl"),
                        seed = TRUE)) %dofuture% {
              
              # To avoid issues with native libraries, try keeping native libraries thread safe
              data.table::setDTthreads(1L)
              RhpcBLASctl::blas_set_num_threads(1L)
              RhpcBLASctl::omp_set_num_threads(1L)
              
              # Generate and analyse data
               data <- cohort_time_var_past_e_u(n = 100000, obs_time = 500, risk_window = risk_win,
                                               p_U = 1e-3, U_on_C = 1, U_on_Y = 1,     
                                               p_C = p_C, C_on_C = 0.3, E_on_C = 0.2,
                                               p_E = 3.2e-3, C_on_E = C_on_E, E_on_E = 1e-7,
                                              baseline_Y = 2e-5, IRR_E = IRR_E, IRR_C = IRR_C
                                                                     )
                                              data_SCCS <- SCCS_reformat_confound(data, risk_window = risk_win)
                                              analyse_sccs_confound(data= data_SCCS, rep = i, C_adjust = FALSE)
                        }
      
  # Export results
      out_file <- file.path(out_directory, paste0(scen_name, ".csv"))
      data.table::fwrite(sim_res, out_file)
                         
  # Evaluate bias and coverage
  perf <- bias_quantification(true_IRR_E = IRR_E, true_IRR_C = IRR_C,
                              result_table = sim_res)
                         
  cbind(scenarios[s, ], perf)
                       }
out_file2 <- file.path(out_directory, paste0("bias_scenario_",scen_start,"-",scen_stop, ".csv"))
data.table::fwrite(results_all, out_file2)
View(results_all)


# Simulation approach 2: Adjusting for C in SCCS model ---------------------------
#################################################################################

## Parameter grid ------------------
p_U_vals    <- c(1e-3, 5e-3)
U_on_Y_vals <- c(0.25, 2, 5, 15)
E_on_C_vals <- c(0.2, 2)
IRR_E_vals  <- c(1, 2, 5)
chosen_risk_window <- c(10, 50)

scenarios_adjC <- expand.grid(p_U = p_U_vals,
                              UonY = U_on_Y_vals,
                              EonC = E_on_C_vals,
                              IRR_E = IRR_E_vals,
                              risk_win = chosen_risk_window)

scenarios_adjC$scenario_id <- seq(1:nrow(scenarios_adjC))
scenarios_adjC$scenario_name <- with(scenarios_adjC, 
                                     paste0("Scen_",scenario_id,"_pU", p_U, "_UonY", UonY,
                                            "_EonC", EonC, "_IRRE",
                                            IRR_E, "_",risk_win, "d"))

data.table::fwrite(scenarios_adjC, here("SCCS_Simulation", "Scenarios_Adjust_C.csv"))

# Set the directory to export files
out_directory <- here::here("SCCS_Simulation", "Results", "Adjust_C2")

## Run the simulation ---------------------------------------------------------

### Scenarios 17 to 32 ---------------------------------------------------------
n_sim <- 1000

plan(multisession, workers = 30)

registerDoRNG()
set.seed(1010)

scen_start <- 17
scen_stop <- 32
tic("Sim: bias when adjusting for C") #16:05 Aug 21
results_all_adjC <- foreach(s = scen_start:scen_stop,
                       .combine = rbind) %do% {
   tryCatch({                         
   p_U   <- scenarios_adjC$p_U[s]
   UonY <- scenarios_adjC$UonY[s]
   EonC <- scenarios_adjC$EonC[s]
   IRR_E <- scenarios_adjC$IRR_E[s]
   risk_win <- scenarios_adjC$risk_win[s]
   scen_name <- scenarios_adjC$scenario_name[s]
                         
      # Run simulation for one scenario
      sim_res <- foreach(i = 1:n_sim,
                        .combine = rbind,
                        .options.future = list(packages = c("extraDistr","dplyr","SCCS","data.table","RhpcBLASctl"),
                        seed = TRUE)) %dofuture% {
              

              # To avoid issues with native libraries, try keeping native libraries thread safe
              data.table::setDTthreads(1L)
              RhpcBLASctl::blas_set_num_threads(1L)
              RhpcBLASctl::omp_set_num_threads(1L)
              
              # Simulate data                                                       
              data <- cohort_time_var_past_e_u(n = 100000, obs_time = 500, risk_window = risk_win,
                                                p_U = p_U, U_on_C = 5, U_on_Y = UonY,     
                                                p_C = 1e-3, C_on_C = 0.3, E_on_C = EonC,
                                                p_E = 3.2e-3, C_on_E = 0.25, E_on_E = 1e-7,
                                                baseline_Y = 2e-5, IRR_E = IRR_E, IRR_C = 2
                                                                     )
              data_SCCS <- SCCS_reformat_confound(data, risk_window = risk_win)
              analyse_sccs_confound(data= data_SCCS, rep = i, C_adjust = TRUE)
                } 
                         
     # Export results
     out_file <- file.path(out_directory, paste0(scen_name, ".csv"))
     data.table::fwrite(sim_res, out_file)

     # Evaluate bias and coverage
     perf <- bias_quantification(true_IRR_E = IRR_E, true_IRR_C = 2,
                                  result_table = sim_res)
                         
     cbind(scenarios_adjC[s, ], perf) 
     }, error = function(e) {
       message(sprintf("Error in iteration %d: %s", i, e$message))
       return(NULL)
     })
   }

out_file2 <- file.path(out_directory, paste0("bias_adjC_scenario_",scen_start,"-",scen_stop, ".csv"))
data.table::fwrite(results_all_adjC, out_file2)
View(results_all_adjC)

toc()

### Scenarios 33 to 48 ------------------------------------------------------
n_sim <- 1000

plan(multisession, workers = 30)

registerDoRNG()
set.seed(1010)

scen_start <- 33
scen_stop <- 48
tic("Sim: bias when adjusting for C") #19:58 Aug 21
results_all_adjC <- foreach(s = scen_start:scen_stop,
                            .combine = rbind) %do% {
                              tryCatch({                         
                                p_U   <- scenarios_adjC$p_U[s]
                                UonY <- scenarios_adjC$UonY[s]
                                EonC <- scenarios_adjC$EonC[s]
                                IRR_E <- scenarios_adjC$IRR_E[s]
                                risk_win <- scenarios_adjC$risk_win[s]
                                scen_name <- scenarios_adjC$scenario_name[s]
                                
                                # Run simulation for one scenario
                                sim_res <- foreach(i = 1:n_sim,
                                                   .combine = rbind,
                                                   .options.future = list(packages = c("extraDistr","dplyr","SCCS","data.table","RhpcBLASctl"),
                                                                          seed = TRUE)) %dofuture% {
                                                                            
                                                                            
                                                                            # To avoid issues with native libraries, try keeping native libraries thread safe
                                                                            data.table::setDTthreads(1L)
                                                                            RhpcBLASctl::blas_set_num_threads(1L)
                                                                            RhpcBLASctl::omp_set_num_threads(1L)
                                                                            
                                                                            # Simulate data                                                       
                                                                            data <- cohort_time_var_past_e_u(n = 100000, obs_time = 500, risk_window = risk_win,
                                                                                                             p_U = p_U, U_on_C = 5, U_on_Y = UonY,     
                                                                                                             p_C = 1e-3, C_on_C = 0.3, E_on_C = EonC,
                                                                                                             p_E = 3.2e-3, C_on_E = 0.25, E_on_E = 1e-7,
                                                                                                             baseline_Y = 2e-5, IRR_E = IRR_E, IRR_C = 2
                                                                            )
                                                                            data_SCCS <- SCCS_reformat_confound(data, risk_window = risk_win)
                                                                            analyse_sccs_confound(data= data_SCCS, rep = i, C_adjust = TRUE)
                                                                          } 
                                
                                # Export results
                                out_file <- file.path(out_directory, paste0(scen_name, ".csv"))
                                data.table::fwrite(sim_res, out_file)
                                
                                # Evaluate bias and coverage
                                perf <- bias_quantification(true_IRR_E = IRR_E, true_IRR_C = 2,
                                                            result_table = sim_res)
                                
                                cbind(scenarios_adjC[s, ], perf) 
                              }, error = function(e) {
                                message(sprintf("Error in iteration %d: %s", i, e$message))
                                return(NULL)
                              })
                            }

out_file2 <- file.path(out_directory, paste0("bias_adjC_scenario_",scen_start,"-",scen_stop, ".csv"))
data.table::fwrite(results_all_adjC, out_file2)
View(results_all_adjC)

toc() # 3.5 hours