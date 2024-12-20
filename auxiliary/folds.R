################################################################################
##### Create folds for training

data_folds <- function(species, k = 10) {
  
  selection <- data.table(species = species,
                    row = 1:length(species))
  
  unique_species <- unique(selection$species)
  
  index <- list()
  
  sp_fold <- function(X, unique_species, selection) {
    
    collector <- as.integer()
    
    for(ii in 1:length(unique_species)) {
      
      sub <- subset(selection, species == unique_species[ii])
      samp <- sample(1:nrow(sub), 1)
      collector <- c(collector, sub$row[-samp])
      
    }
    
    return(collector)
    
  }
  
  folds <- lapply(X = 1:k, 
                  FUN = sp_fold, 
                  unique_species= unique_species,
                  selection = selection)
  
  names(folds) <- paste0("Fold", gsub(" ", "0", format(seq(along = 1:k))))
  
  return(folds)
  
}
