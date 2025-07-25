#-------------------------------------------------------------------------------
# Function for build the final models based on the optimal number of components

# PLSDA
model_build_plsda <- function(meta, 
                              split, 
                              segments,
                              rank,
                              spectra,
                              ncomp,
                              threads) {
  
  # Capture the column name
  rank_name <- deparse(substitute(rank))
  
  # Check column exists
  if (!rank_name %in% names(meta)) {
    stop(paste("Column", rank_name, "not found in metadata."))
  }
  
  # Extract the taxon label vector
  species <- meta[[rank_name]]
  
  # Data for training
  meta_split <- meta[split, ]
  species_split <- species[split]
  spectra_split <- spectra[split, ]
  
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"
  
  # Application loop
  build <- function(X, 
                    split, 
                    segments,
                    meta_split,
                    frame,
                    ncomp) {
    
    segment_idx <- segments[[X]]
    relative_idx <- which(split %in% segment_idx)
    
    if (length(relative_idx) < 2) {
      return(try(stop("Too few samples in segment"), silent = TRUE))
    }
    
    meta_segments <- meta_split[relative_idx, ]
    sub_frame <- frame[relative_idx, ]
    
    if (length(unique(sub_frame$species)) < 2) {
      return(try(stop("Only one class in segment"), silent = TRUE))
    }
    
    ctrl <- trainControl(index = data_folds(sub_frame$species, k = 10))
    
    model <- try(train(species ~ .,
                       data = sub_frame,
                       method = "pls",
                       trControl = ctrl,
                       preProcess = NULL,
                       tuneGrid = expand.grid(ncomp = ncomp),
                       metric = "Accuracy",
                       maxit = 10000),
                 silent = TRUE)
    
    return(model)
  }
  
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
                            rank,
                            spectra,
                            threads) {
  
  # Resolve the column name from unquoted argument
  rank_name <- deparse(substitute(rank))
  
  # Check that column exists
  if (!rank_name %in% names(meta)) {
    stop(paste("Column", rank_name, "not found in metadata."))
  }
  
  # Extract species/taxon vector
  species <- meta[[rank_name]]
  
  # Subset to split
  meta_split <- meta[split, ]
  species_split <- species[split]
  spectra_split <- spectra[split, ]
  
  # Build data frame for modeling
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"
  
  # One model per segment
  build <- function(X, split, segments, meta_split, frame) {
    segment_idx <- segments[[X]]
    
    # Map global sample indices to local split positions
    relative_idx <- which(split %in% segment_idx)
    
    if (length(relative_idx) < 2) {
      return(try(stop("Too few samples in segment"), silent = TRUE))
    }
    
    sub_frame <- frame[relative_idx, ]
    
    # Corrected: check class diversity using sub_frame, not meta_segment
    if (length(unique(sub_frame$species)) < 2) {
      return(try(stop("Only one class in segment"), silent = TRUE))
    }
    
    # Generate CV folds (you can replace data_folds if needed)
    ctrl <- trainControl(index = data_folds(sub_frame$species, k = 10))
    
    model <- try(train(species ~ .,
                       data = sub_frame,
                       method = "lda",
                       trControl = ctrl,
                       preProcess = NULL,
                       metric = "Accuracy"),
                 silent = TRUE)
    
    return(model)
  }
  
  # Run in parallel
  complete <- pbmclapply(X = seq_along(segments), 
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
