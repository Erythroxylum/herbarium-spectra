#-------------------------------------------------------------------------------
# Function for data split based on species and accession

data_split <- function(meta) {
  
  unique_species <- unique(meta$scientificName)
  selected_rows <- integer()
  
  for (i in seq_along(unique_species)) {
    
    # Subset metadata by species
    sub_meta <- subset(meta, scientificName == unique_species[i])
    
    # Get all specimen identifiers for that species
    unique_specimens <- unique(sub_meta$specimenIdentifier)
    
    # Randomly select up to 10 specimens (adjust if fewer available)
    selected_ids <- sample(unique_specimens, min(10, length(unique_specimens)))
    
    # Subset rows corresponding to selected specimens
    sub_rows <- which(meta$specimenIdentifier %in% selected_ids & meta$scientificName == unique_species[i])
    
    # Collect row indices
    selected_rows <- c(selected_rows, sub_rows)
  }
  
  return(selected_rows)
}

