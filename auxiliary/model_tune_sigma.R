#-------------------------------------------------------------------------------
# Function for tune the optimal number of components of the PLSR models

model_tune <- function(meta,  
                       split,  
                       segments,
                       traits,  
                       spectra,
                       ncomp_max,
                       threads) {  
  # Data for training
  meta_split <- meta[split,]
  traits_split <- traits[split,]
  spectra_split <- spectra[split,]
  
  frame <- cbind(traits_split, spectra_split)
  colnames(frame)[1] <- "trait"
  
  # Application loop
  tune <- function(X,  
                   split,  
                   segments,
                   meta_split,
                   frame,
                   ncomp_max) {  
    # Get samples from segments to train
    get_segments <- split %in% segments[[X]]
    meta_segments <- meta_split[get_segments,]
    sub_frame <- frame[get_segments,]
    #kfolds <- data_folds(meta_segments$species, k = 10)
    
    # Perform model
    ml_model <- plsr(trait ~ .,  
                     data = sub_frame,
                     validation = "CV",
                     ncomp = ncomp_max,
                     center = TRUE,
                     method = "oscorespls",
                     maxit = 10000,
                     segments = 10) # was kfolds for balanced sampling
    
    # Calculate optimal number of components using selectNcomp
    ncomp <- selectNcomp(ml_model, method = "onesigma", plot = FALSE)
    
    # Performance
    RMSEP_train <- pls::RMSEP(object = ml_model,
                              estimate = "train")$val[,,]
    
    RMSEP_CV <- pls::RMSEP(object = ml_model,
                           estimate = "CV")$val[,,]
    
    R2_train <- pls::R2(object = ml_model,
                        estimate = "train")$val[,,]
    
    R2_CV <- pls::R2(object = ml_model,
                     estimate = "CV")$val[,,]
    
    PRESS <- c(Intercept = 0, ml_model$validation$PRESS[1,])
    
    # Collect results
    results <- data.table(iteration = X,
                          nsamples = nrow(sub_frame),
                          ncomp = ncomp,  # Add ncomp to the results
                          metric = c("RMSEP", "RMSEP", "R2", "R2", "PRESS"),
                          estimate = c("train", "CV", "train", "CV", "PRESS"))
    metrics <- rbind(RMSEP_train, RMSEP_CV, R2_train, R2_CV, PRESS)
    results <- cbind(results, metrics)
    
    # Clean memory
    rm(list = c("get_segments", "sub_frame", "ml_model",  
                "RMSEP_train", "RMSEP_CV", "R2_train", "R2_CV",  
                "PRESS", "metrics"))
    
    return(results)
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:length(segments),  
                         FUN = tune,  
                         meta_split = meta_split,
                         split = split,  
                         segments = segments,
                         frame = frame,
                         ncomp_max = ncomp_max,
                         mc.preschedule = TRUE,  
                         mc.set.seed = FALSE,
                         mc.cores = threads,
                         mc.cleanup = TRUE,  
                         mc.allow.recursive = TRUE)
  
  frame <- data.table()
  
  for (i in 1:length(complete)) {
    frame <- rbind(frame, complete[[i]], fill = TRUE)
  }
  
  return(frame)
}
