###################################
### Helpder functions for       ###
### creating directories and    ###
### saving results.             ###
###################################

# Create the folders to store results
create_directory <- function(path) {
  if (!dir.exists(path)) 
    dir.create(path, recursive = TRUE)
}

# Saving results to .csv file 
append_to_csv <- function(out_df, file_path) {
  if (!file.exists(file_path)) {
    write.csv(out_df, file_path, row.names = FALSE)
  } else {
    write.table(out_df, file_path, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE)
  }
}

