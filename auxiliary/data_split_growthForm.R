#-------------------------------------------------------------------------------
# Function for data split based on growth form following Kothari et al 2022

# Function for data split based on growth form and accession
data_split <- function(meta, p = 0.75) {
  
  unique_growthforms <- unique(meta$growthForm)
  collector <- numeric()
  
  for (gf in unique_growthforms) {
    
    # Subset metadata for the current GrowthForm
    sub_meta <- subset(meta, growthForm == gf)
    
    # Determine unique accessions and the number to include in the training set
    unique_accessions <- unique(sub_meta$accession)
    n_accessions <- floor(p * length(unique_accessions))
    
    # Randomly select accessions for the training set
    selected_accessions <- sample(unique_accessions, n_accessions)
    
    # Include all indices corresponding to the selected accessions
    train_indices <- which(meta$accession %in% selected_accessions)
    
    # Save the indices to the collector
    collector <- c(collector, train_indices)
  }
  
  return(collector)
}

#### old

data_split_allIndices <- function(meta, p = 0.75) {
  
  unique_growthforms <- unique(meta$growthForm)
  collector <- numeric()
  
  for (gf in unique_growthforms) {
    
    # Subset metadata for the current GrowthForm
    sub_meta <- subset(meta, growthForm == gf)
    
    # Determine the number of samples to include in the training set
    n_samples <- floor(p * nrow(sub_meta))
    
    # Randomly select indices for the training set
    train_indices <- sample(seq_len(nrow(sub_meta)), n_samples)
    
    # Save the selected training samples to the collector
    collector <- c(collector, sub_meta$sample[train_indices])
  }
  
  return(collector)
}

