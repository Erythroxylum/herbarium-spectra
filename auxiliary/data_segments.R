#-------------------------------------------------------------------------------
# Function for creating segments of data for training based on accessions and leaves

data_segments <- function(X, 
                          meta,
                          split) {
  
  # Select samples on training
  meta_split <- meta[split,]
  
  # Unique accession
  unique_accession <- unique(meta_split$specimenIdentifier)
  
  # Collector
  collector <- as.numeric()
  
  # Loop for selection
  for(i in 1:length(unique_accession)) {
    
    sub_accession <- subset(meta_split, specimenIdentifier == unique_accession[i])
    samp <- sample(sub_accession$sample, 1)
    collector <- c(collector, samp)
    
  }
  
  return(collector)
  
}
