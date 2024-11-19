#-------------------------------------------------------------------------------
# Function for data split based on species and accession

data_split <- function(meta) {
  
  unique_species <- unique(meta$species)
  collector <- as.numeric()
  
  for(i in 1:length(unique_species)) {
    
    # Subset species of interest and look to the accession
    sub_meta <- subset(meta, species == unique_species[i])
    unique_accession <- unique(sub_meta$accession)
    
    # Ramdomly select accession
    get_samples <- sample(1:length(unique_accession), 10)
    
    # Select accession to subset for traint
    get_samples <- sub_meta$accession %in% unique_accession[get_samples]
    sub_meta <- sub_meta[get_samples == TRUE, ]
    
    # Save collector of sample
    collector <- c(collector, sub_meta$sample)
    
  }
  
  return(collector)
  
}

