#-------------------------------------------------------------------------------
# Function for tune the optimal number of components of the PLSR models

model_tune_plsda <- function(meta,
                       split,
                       segments,
                       rank,
                       spectra,
                       ncomp_max,
                       threads) {

  # Data for training
  meta_split <- meta[split,]
  traits_split <- traits[split,]
  rank_split <- rank[split,] # [split,5] for testing
  spectra_split <- spectra[split,]

  # combine taxon meta with spectra, then rename the taxon level to test as "taxon"
  #frame <- cbind(traits_split, spectra_split)
  frame <- cbind(rank_split, spectra_split)
  ###### should the other levels of rank be removed?

  # Application loop
  tune <- function(X,
                   split,
                   segments,
                   frame,
                   ncomp_max=30) {

    # Get samples from segments to train
    get_segments <- split %in% segments[[X]]
    sub_frame <- frame[get_segments,]

    # Perform model
    #ml_model <- plsr(rank ~ .,
     #                data = sub_frame,
      #               validation = "CV",
       #              ncomp = ncomp_max,
        #             center = TRUE,
         #            method = "oscorespls",
          #           maxit = 10000)

    ctrl <- trainControl(method = "cv",
                         #repeats = 10,
                         number = 10,
                         summaryFunction = multiClassSummary,
                         #sampling = "none"
                         classProbs = TRUE
                         )

    pls_model <- train(spShort ~ .,
                       data = sub_frame,
                       method = "pls",
                       trControl = ctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 30,
                       metric = 'ROC')

    # Performance
    predictions <- predict(pls_model, newdata = as.matrix(sub_frame), type = "prob")

    RMSEP_train <- postResample(predict(pls_model$bestTune, as.matrix(spectra_split)), sub_frame$spShort)
    RMSEP_CV <- ml_model$resample$RMSE
    R2_train <- cor(predict(ml_model, sub_frame), sub_frame$rank)^2  # R^2 necessary?

    #RMSEP_train <- pls::RMSEP(object = ml_model,
     #                         estimate = "train")$val[,,]

    #RMSEP_CV <- pls::RMSEP(object = ml_model,
     #                      estimate = "CV")$val[,,]

    #R2_train <- pls::R2(object = ml_model,
     #                   estimate = "train")$val[,,]

    #R2_CV <- pls::R2(object = ml_model,
     #                estimate = "CV")$val[,,]

    #PRESS <- c(Intercept = 0, ml_model$validation$PRESS[1,])

    # Collect results
    results <- data.table(iteration = X,
                          nsamples = nrow(sub_frame),
                          metric = c("RMSEP", "R2"),
                          estimate = c("train", "train"))

    #results <- data.table(iteration = X,
     #                     nsamples = nrow(sub_frame),
      #                    metric = c("RMSEP", "RMSEP", "R2", "R2", "PRESS"),
       #                   estimate = c("train", "CV", "train", "CV", "PRESS"))


    metrics <- data.table(RMSEP = RMSEP_train[1], R2 = R2_train)
    results <- cbind(results, metrics)

    #metrics <- rbind(RMSEP_train, RMSEP_CV, R2_train, R2_CV, PRESS)
    #results <- cbind(results, metrics)

    # Clean memory
    rm(list = c("get_segments", "sub_frame", "ml_model", "metrics"))
    #rm(list = c("get_segments", "sub_frame", "ml_model",
     #           "RMSEP_train", "RMSEP_CV", "R2_train", "R2_CV",
      #          "PRESS", "metrics"))

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
