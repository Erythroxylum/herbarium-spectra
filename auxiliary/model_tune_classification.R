#-------------------------------------------------------------------------------
# Function for tune the optimal number of components of the PLSR models

model_tune_plsda <- function(meta,
                             split,
                             segments,
                             species,
                             spectra,
                             ncomp_max,
                             threads) {


  # Data for training
  meta_split <- meta[split,]
  species_split <- as.factor(species[split])
  spectra_split <- spectra[split,]
  
  # Combine taxa with spectra
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"


  # Application loop
  tune <- function(X,
                   split,
                   segments,
                   frame,
                   meta,
                   ncomp_max) {

    # Get samples from segments to train
    get_segments <- split %in% segments[[X]]
    meta_split <- meta_split[get_segments,]
    sub_frame <- frame[get_segments,]

    # Model control
    ctrl <- trainControl(index = data_folds(meta_split$species, k = 10))

    model <- train(species ~ .,
                   data = sub_frame,
                   method = "pls",
                   trControl = ctrl,
                   preProcess = c("center"),
                   tuneGrid = expand.grid(ncomp = 1:ncomp_max),
                   metric = "Accuracy",
                   maxit = 10000)

    # Performance
    results <- data.table(iteration = X,
                          nsamples = nrow(sub_frame),
                          metric = colnames(model$results)[-1])

    results <- cbind(results, t(model$results)[-1,])
    
    # Clean memory
    rm(list = c("get_segments", "sub_frame", "model"))

    return(results)

  }

  # Parallel processing
  complete <- pbmclapply(X = 1:length(segments),
                         FUN = tune,
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

  for(i in 1:length(complete)) {
    frame <- rbind(frame, complete[[i]], fill = TRUE)
  }

  return(frame)

}
