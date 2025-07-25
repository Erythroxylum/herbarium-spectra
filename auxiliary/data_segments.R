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

## For classification

data_segments_classification <- function(X, meta, split) {
  repeat {
    # Filter training data
    meta_split <- meta[split, ]
    
    # Sample one measurement per specimen
    sampled_rows <- meta_split[, .SD[sample(.N, 1)], by = specimenIdentifier]$sample
    
    # Check that there is class diversity
    if (length(unique(meta[sampled_rows]$scientificName)) >= 2) {
      return(sampled_rows)
    }
    # Else: repeat until segment contains at least 2 species
  }
}
