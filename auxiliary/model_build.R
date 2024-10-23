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
                    frame,
                    ncomp) {
    
    # Get samples from segments to train
    get_segments <- split %in% segments[[X]]
    sub_frame <- frame[get_segments,]
    
    # Perform model
    ml_model <- plsr(trait ~ ., 
                     data = sub_frame,
                     validation = "CV",
                     ncomp = ncomp,
                     center = TRUE,
                     method = "oscorespls",
                     maxit = 10000)
    
    
    
    return(ml_model)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:length(segments), 
                         FUN = build, 
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
