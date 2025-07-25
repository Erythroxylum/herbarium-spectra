#-------------------------------------------------------------------------------
# Function for tune the optimal number of components of the PLSR models

model_tune_plsda <- function(meta,
                             split,
                             segments,
                             rank,
                             spectra,
                             ncomp_max,
                             threads) {
  
  # Extract species/taxon name from rank argument
  rank_name <- deparse(substitute(rank))
  if (!rank_name %in% names(meta)) {
    stop(paste("Column", rank_name, "not found in metadata."))
  }
  
  species <- meta[[rank_name]]
  
  # Data for training
  meta_split <- meta[split, ]
  species_split <- factor(species[split])
  spectra_split <- spectra[split, ]
  
  # Combine taxa with spectra
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"
  
  # Application loop
  tune <- function(X, split, segments, frame, meta, ncomp_max) {
    segment_idx <- segments[[X]]
    relative_idx <- which(split %in% segment_idx)
    
    if (length(relative_idx) < 2) {
      return(NULL)
    }
    
    meta_subset <- meta[split, ][relative_idx, ]
    sub_frame <- frame[relative_idx, ]
    
    # Use class labels from the modeling frame
    ctrl <- trainControl(index = data_folds(sub_frame$species, k = 10))
    
    model <- train(species ~ .,
                   data = sub_frame,
                   method = "pls",
                   trControl = ctrl,
                   preProcess = NULL,
                   tuneGrid = expand.grid(ncomp = 1:ncomp_max),
                   metric = "Accuracy",
                   maxit = 10000)
    
    # Format results
    results <- data.table(
      iteration = X,
      nsamples = nrow(sub_frame),
      metric = colnames(model$results)[-1]
    )
    results <- cbind(results, t(model$results)[-1, ])
    
    return(results)
  }
  
  # Run in parallel
  complete <- pbmclapply(
    X = seq_along(segments),
    FUN = tune,
    split = split,
    segments = segments,
    frame = frame,
    meta = meta,
    ncomp_max = ncomp_max,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE,
    mc.cores = threads,
    mc.cleanup = TRUE,
    mc.allow.recursive = TRUE
  )
  
  # Combine into one table
  frame_out <- rbindlist(complete, fill = TRUE)
  return(frame_out)
}
