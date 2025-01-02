#-------------------------------------------------------------------------------
# Function for build the final models based on the optimal number of components

# PLSDA
model_build_plsda <- function(meta, 
                              split, 
                              segments,
                              species, 
                              spectra,
                              ncomp,
                              threads) {
  
  # Data for training
  meta_split <- meta[split,]
  species_split <- species[split]
  spectra_split <- spectra[split,]
  
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"
  
  # Application loop
  build <- function(X, 
                    split, 
                    segments,
                    meta_split,
                    frame,
                    ncomp) {
    
    # Get samples from segments to train
    get_segments <- split %in% segments[[X]]
    meta_split <- meta_split[get_segments,]
    sub_frame <- frame[get_segments,]
    
    # Model control
    ctrl <- trainControl(index = data_folds(meta_split$species, k = 10))
    
    # Perform model
    model <- train(species ~ .,
                   data = sub_frame,
                   method = "pls",
                   trControl = ctrl,
                   preProcess = NULL,
                   tuneGrid = expand.grid(ncomp = ncomp),
                   metric = "Accuracy",
                   maxit = 10000)
    
    # Clean memory
    rm(list = c("get_segments", "meta_split", "sub_frame"))
    
    return(model)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:length(segments), 
                         FUN = build, 
                         split = split, 
                         segments = segments,
                         meta_split = meta_split,
                         frame = frame,
                         ncomp = ncomp,
                         mc.preschedule = TRUE, 
                         mc.set.seed = FALSE,
                         mc.cores = threads,
                         mc.cleanup = TRUE, 
                         mc.allow.recursive = TRUE)
  
  return(complete)
  
}

# LDA
model_build_lda <- function(meta, 
                            split, 
                            segments,
                            species, 
                            spectra,
                            threads) {
  
  # Data for training
  meta_split <- meta[split,]
  species_split <- species[split]
  spectra_split <- spectra[split,]
  
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"
  
  # Application loop
  build <- function(X, 
                    split, 
                    segments,
                    meta_split,
                    frame) {
    
    # Get samples from segments to train
    get_segments <- split %in% segments[[X]]
    meta_split <- meta_split[get_segments,]
    sub_frame <- frame[get_segments,]
    
    # Model control
    ctrl <- trainControl(index = data_folds(meta_split$species, k = 10))
    
    # Perform model
    model <- train(species ~ .,
                   data = sub_frame,
                   method = "lda",
                   trControl = ctrl,
                   preProcess = NULL,
                   metric = "Accuracy",
                   maxit = 10000)
    
    # Clean memory
    rm(list = c("get_segments", "meta_split", "sub_frame"))
    
    return(model)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:length(segments), 
                         FUN = build, 
                         split = split, 
                         segments = segments,
                         meta_split = meta_split,
                         frame = frame,
                         mc.preschedule = TRUE, 
                         mc.set.seed = FALSE,
                         mc.cores = threads,
                         mc.cleanup = TRUE, 
                         mc.allow.recursive = TRUE)
  
  return(complete)
  
}
