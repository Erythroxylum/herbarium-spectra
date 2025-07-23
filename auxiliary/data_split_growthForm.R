#-------------------------------------------------------------------------------
# Function for data split based on growth form following Kothari et al 2022

# Function for data split based on growth form and specimenIdentifier
data_split <- function(meta, p = 0.75) {
  
  unique_growthforms <- unique(meta$growthForm)
  collector <- numeric()
  
  for (gf in unique_growthforms) {
    
    # Subset metadata for the current GrowthForm
    sub_meta <- subset(meta, growthForm == gf)
    
    # Determine unique specimenIdentifiers and the number to include in the training set
    unique_specimenIdentifiers <- unique(sub_meta$specimenIdentifier)
    n_specimenIdentifiers <- floor(p * length(unique_specimenIdentifiers))
    
    # Randomly select specimenIdentifiers for the training set
    selected_specimenIdentifiers <- sample(unique_specimenIdentifiers, n_specimenIdentifiers)
    
    # Include all indices corresponding to the selected specimenIdentifiers
    train_indices <- which(meta$specimenIdentifier %in% selected_specimenIdentifiers)
    
    # Save the indices to the collector
    collector <- c(collector, train_indices)
  }
  
  return(collector)
}
