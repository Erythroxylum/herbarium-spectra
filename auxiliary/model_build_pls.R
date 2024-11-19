#-------------------------------------------------------------------------------
# Function for build the final models based on the optimal number of components

model_build <- function(meta, 
                        split, 
                        segments,
                        traits, 
                        spectra,
                        ncomp,
                        threads) {
  
  # Data for training
  meta_split <- meta[split,]
  traits_split <- traits[split,]
  spectra_split <- spectra[split,]
  
  frame <- cbind(traits_split, spectra_split)
  colnames(frame)[1] <- "trait"
  
  # Application loop
  build <- function(X, 
                    split, 
                    segments,
                    meta_split,
                    frame,
                    ncomp) {
    
    # Get samples from segments to train
    get_segments <- split %in% segments[[X]]
    meta_segments <- meta_split[get_segments,]
    sub_frame <- frame[get_segments,]
    #kfolds <- data_folds(meta_segments$species, k = 10)
    
    # Perform model
    ml_model <- plsr(trait ~ ., 
                     data = sub_frame,
                     validation = "CV",
                     ncomp = ncomp,
                     center = TRUE,
                     method = "oscorespls",
                     maxit = 10000,
                     segments = 10) # was kfolds for balanced sampling
    
    return(ml_model)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:length(segments), 
                         FUN = build, 
                         meta_split = meta_split,
                         split = split, 
                         segments = segments,
                         frame = frame,
                         ncomp = ncomp,
                         mc.preschedule = TRUE, 
                         mc.set.seed = FALSE,
                         mc.cores = threads,
                         mc.cleanup = TRUE, 
                         mc.allow.recursive = TRUE)
  
  return(complete)
  
}
