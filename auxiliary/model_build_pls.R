#-------------------------------------------------------------------------------
# Function for build the final models based on the optimal number of components

model_build <- function(meta,
                        split,
                        segments,
                        traits,
                        spectra,
                        ncomp,
                        threads = 1) {
  
  # Ensure proper types
  traits_split <- traits[split][[1]]
  spectra_split <- as.data.frame(spectra[split, ])
  meta_split <- meta[split, ]
  frame <- as.data.frame(cbind(trait = traits_split, spectra_split))
  
  # Inner model-building function
  build <- function(X, split, segments, meta_split, frame, ncomp) {
    get_segments <- split %in% segments[[X]]
    sub_frame <- frame[get_segments, ]
    meta_seg <- meta_split[get_segments, ]
    
    if (nrow(sub_frame) < ncomp) {
      stop(paste0("Segment ", X, ": not enough samples to fit ", ncomp, " components"))
    }
    
    plsr(trait ~ .,
         data = sub_frame,
         validation = "CV",
         ncomp = ncomp,
         center = TRUE,
         method = "oscorespls",
         maxit = 10000,
         segments = 10)
  }
  
  # Set up parallel plan and progress bar
  plan(multisession, workers = threads)
  handlers(global = TRUE)
  handlers("txtprogressbar")  # You can change this to "progress" or "notifier"
  
  # Wrap parallel processing in progressr
  p <- progressr::progressor(along = 1:length(segments))
  
  complete <- future_lapply(1:length(segments), function(i) {
    p(sprintf("Fitting model %d", i))
    build(i, split, segments, meta_split, frame, ncomp)
  }, future.seed = TRUE)
  
  return(complete)
}

# OLD
model_build_old <- function(meta, 
                        split, 
                        segments,
                        traits, 
                        spectra,
                        ncomp,
                        threads = 1) {
  
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
