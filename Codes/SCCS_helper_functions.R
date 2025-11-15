##########################################
### Project: Time-varying confounding  ###
### affected by past exposure in SCCS  ###
##########################################

###################################
### Helper functions for        ###
### running the simulation      ###
###################################

# Create the folders to store results ---------------------------------------
create_directory <- function(path) {
  if (!dir.exists(path)) 
    dir.create(path, recursive = TRUE)
}

# Saving results to .csv file -----------------------------------------------
append_to_csv <- function(out_df, file_path) {
  if (!file.exists(file_path)) {
    write.csv(out_df, file_path, row.names = FALSE)
  } else {
    write.table(out_df, file_path, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE)
  }
}

# Log errors ----------------------------------------------------------------

log_error <- function(e, stage, seed = NA, 
                      scen_name = NA, rep = NA,
                      method = NA, 
                      output_dir = here("Results")) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  log_file <- file.path(output_dir, "error_log.txt")
  
  # Log message
  msg <- paste0(
    "[", Sys.time(), "] ",
    "Stage: ", stage, " | ",
    "Scenario: ", scen_name, " | ",
    "Rep: ", rep, " | ",
    "Method: ", method, " | ",
    "Seed: ", seed, "\n",
    "Error message: ", conditionMessage(e), "\n\n"
  )
  
  # Append message to log file
  cat(msg, file = log_file, append = TRUE)
}



