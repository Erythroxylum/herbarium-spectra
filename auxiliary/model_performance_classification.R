#-------------------------------------------------------------------------------
# Model Performance Function for PLS-DA
model_performance_plsda <- function(meta_split,
                                    species_split, 
                                    spectra_split,
                                    models,
                                    ncomp,
                                    threads = 1) {
  
  # Clean species labels and prepare prediction frame
  species_split <- gsub("_", " ", species_split)
  frame <- cbind(species = species_split, spectra_split)
  
  #-----------------------------------------------------------------------------
  # Predict function (safe)
  application_predict <- function(X, models, frame, ncomp) {
    predicted <- try(
      predict(models[[X]], newdata = frame, ncomp = ncomp, type = "prob"),
      silent = TRUE
    )
    
    if (inherits(predicted, "try-error") || is.null(predicted)) {
      stop(sprintf("prediction failed for model[%d]", X))
    }
    
    if (!is.matrix(predicted) && !is.data.frame(predicted)) {
      stop(sprintf("prediction for model[%d] is not a matrix or data.frame", X))
    }
    
    predicted <- as.data.table(cbind(iteration = X, predicted))
    return(predicted)
  }
  
  # Run predictions in parallel
  predicted <- pbmclapply(
    X = seq_along(models),
    FUN = application_predict,
    models = models,
    frame = frame,
    ncomp = ncomp,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE,
    mc.cores = threads,
    mc.cleanup = TRUE
  )
  
  #-----------------------------------------------------------------------------
  # Performance evaluation function
  application_performance <- function(X, species_split, prediction) {
    if (is.null(prediction) || !is.data.frame(prediction)) {
      stop(sprintf("prediction[[%d]] is invalid or NULL", X))
    }
    
    prediction_matrix <- prediction[, -1, with = FALSE]  # remove iteration column
    species_names <- colnames(prediction_matrix)
    
    predicted_species <- apply(prediction_matrix, 1, which.max)
    predicted_species <- species_names[predicted_species]
    
    # Ensure consistent factor levels for confusion matrix
    predicted_species <- as.character(predicted_species)
    species_split <- as.character(species_split)
    all_classes <- union(unique(species_split), unique(predicted_species))
    predicted_species <- factor(predicted_species, levels = all_classes)
    species_split <- factor(species_split, levels = all_classes)
    
    # Confusion matrix
    tab <- table(Predicted = predicted_species, Actual = species_split)
    
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      stop(sprintf("Confusion matrix is degenerate at iteration %d", X))
    }
    
    cm <- confusionMatrix(tab)
    
    # Extract performance
    overall <- as.data.table(t(cm$overall))
    byClass <- as.data.table(cm$byClass)
    byClass[, Class := rownames(cm$byClass)]
    stopifnot("Class" %in% names(byClass))
    
    result <- cbind(data.table(iteration = X, ncomp = ncomp), overall, byClass)
    return(result)
  }
  
  # Run performance evaluation in parallel
  performance <- pbmclapply(
    X = seq_along(models),
    FUN = function(X) application_performance(
      X,
      species_split = species_split,
      prediction = predicted[[X]]
    ),
    mc.preschedule = TRUE,
    mc.set.seed = FALSE,
    mc.cores = threads,
    mc.cleanup = TRUE
  )
  
  # Bind results
  performance <- rbindlist(performance, fill = TRUE)
  performance$Class <- gsub("Class: ", "", performance$Class)
  
  predicted <- rbindlist(predicted, fill = TRUE)
  
  return(list(
    performance = performance,
    predicted = predicted,
    ncomp = ncomp
  ))
}




#-------------------------------------------------------------------------------
# Model Performance Function for LDA
model_performance_lda <- function(meta_split,
                                  species_split, 
                                  spectra_split,
                                  models,
                                  threads = 1) {
  
  # Prepare prediction frame
  species_split <- gsub("_", " ", species_split)
  frame <- cbind(species = species_split, spectra_split)
  
  #-----------------------------------------------------------------------------
  # Predict function
  application_predict <- function(X, models, frame) {
    predicted <- try(predict(models[[X]], newdata = frame, type = "prob"), silent = TRUE)
    
    if (inherits(predicted, "try-error") || is.null(predicted)) {
      stop(sprintf("Prediction failed for model[%d]", X))
    }
    
    if (!is.matrix(predicted) && !is.data.frame(predicted)) {
      stop(sprintf("prediction for model[%d] is not a matrix or data.frame", X))
    }
    
    predicted <- as.data.table(cbind(iteration = X, predicted))
    return(predicted)
  }
  
  # Run predictions
  predicted <- pbmclapply(
    X = seq_along(models),
    FUN = application_predict,
    models = models,
    frame = frame,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE,
    mc.cores = threads,
    mc.cleanup = TRUE
  )
  
  #-----------------------------------------------------------------------------
  # Evaluate performance
  application_performance <- function(X, species_split, prediction) {
    if (is.null(prediction) || !is.data.frame(prediction)) {
      stop(sprintf("prediction[[%d]] is invalid or NULL", X))
    }
    
    prediction_matrix <- prediction[, -1, with = FALSE]
    species_names <- colnames(prediction_matrix)
    
    predicted_species <- apply(prediction_matrix, 1, which.max)
    predicted_species <- species_names[predicted_species]
    
    predicted_species <- as.character(predicted_species)
    species_split <- as.character(species_split)
    
    # Match factor levels
    all_classes <- union(unique(species_split), unique(predicted_species))
    predicted_species <- factor(predicted_species, levels = all_classes)
    species_split <- factor(species_split, levels = all_classes)
    
    tab <- table(Predicted = predicted_species, Actual = species_split)
    
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      stop(sprintf("Confusion matrix is degenerate at iteration %d", X))
    }
    
    cm <- confusionMatrix(tab)
    
    overall <- as.data.table(t(cm$overall))
    byClass <- as.data.table(cm$byClass)
    byClass[, Class := rownames(cm$byClass)]
    stopifnot("Class" %in% names(byClass))
    
    result <- cbind(data.table(iteration = X), overall, byClass)
    return(result)
  }
  
  # Evaluate in parallel
  performance <- pbmclapply(
    X = seq_along(models),
    FUN = function(X) application_performance(
      X,
      species_split = species_split,
      prediction = predicted[[X]]
    ),
    mc.preschedule = TRUE,
    mc.set.seed = FALSE,
    mc.cores = threads,
    mc.cleanup = TRUE
  )
  
  performance <- rbindlist(performance, fill = TRUE)
  performance$Class <- gsub("Class: ", "", performance$Class)
  
  predicted <- rbindlist(predicted, fill = TRUE)
  
  return(list(
    performance = performance,
    predicted = predicted
  ))
}
