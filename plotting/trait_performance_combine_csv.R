###########################################################
#### Herbarium or Pressed traits

# Load necessary library
library(dplyr)

# Specify the folder containing the CSV files
folder_path <- paste0(getwd(),"/results/Table_performance_allTraits/")

# Get a list of all CSV files in the folder
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty list to store results
results_list <- list()

# Loop through each file and compute the mean and standard deviation for each column
for (file in file_list) {
  # Read the CSV file
  data <- read.csv(file)
  
  # Compute mean and standard deviation for each column
  stats <- data %>%
    summarise(across(everything(), list(mean = ~mean(.x, na.rm = TRUE),
                                        sd = ~sd(.x, na.rm = TRUE))))
  
  # Flatten the results into a single row
  flat_stats <- unlist(stats)
  
  # Store results in a list, with file name (without extension) as the name
  file_name <- tools::file_path_sans_ext(basename(file))
  results_list[[file_name]] <- flat_stats
}

# Combine all results into a single data frame
combined_results <- do.call(rbind, results_list)

# Convert to data frame for cleaner display
combined_results <- as.data.frame(combined_results)

# Add row names as a column (if needed)
combined_results$file_name <- rownames(combined_results)

# View the final combined results
print(combined_results)

write.csv(combined_results, "results/Table-S2_performance_allTraits.csv")



###########################################################
#### Model transfer

library(dplyr)

# Specify the folder containing the CSV files
folder_path <- paste0(getwd(), "/results/model_transfer/")

# Get a list of all CSV files ending with "*performance.csv" in all subdirectories
file_list <- list.files(path = folder_path, pattern = "performance\\.csv$", full.names = TRUE, recursive = TRUE)

# Initialize an empty list to store results
results_list <- list()

# Loop through each file and compute the mean and standard deviation for each column
for (file in file_list) {
  # Read the CSV file
  data <- read.csv(file)
  
  # Compute mean and standard deviation for each column
  stats <- data %>%
    summarise(across(everything(), list(mean = ~mean(.x, na.rm = TRUE),
                                        sd = ~sd(.x, na.rm = TRUE))))
  
  # Flatten the results into a single row
  flat_stats <- unlist(stats)
  
  # Add the filename to the results
  file_name <- basename(file)
  results_list[[file_name]] <- c(flat_stats, file_name = file_name)
}

# Combine all results into a single data frame
combined_results <- do.call(rbind, results_list)

# Convert to data frame for cleaner display
combined_results <- as.data.frame(combined_results, stringsAsFactors = FALSE)

# View the final combined results
print(combined_results)

# Write the combined results to a CSV file
write.csv(combined_results, "results/LMA_modelTransfer_performance_results.csv", row.names = FALSE)





