#-------------------------------------------------------------------------------
# Model Performance Function for PLSDA
model_performance_plsda <- function(meta_split,
                                    species_split, 
                                    spectra_split,
                                    models,
                                    ncomp,
                                    threads = 1) {
  
  # Data for predicting
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"
  
  #-----------------------------------------------------------------------------
  # Apply prediction function
  application_predict <- function(X,
                                  models,
                                  frame,
                                  ncomp) {
    # Predict new values
    predicted <- predict(models[[X]], 
                         newdata = frame,
                         ncom = ncomp,
                         type = "prob")
    
    predicted <- as.data.table(cbind(iteration = X,
                                     predicted))
    
    #predicted <- as.data.table(predicted)
    
    return(predicted)
    
  }
  
  # Predicted in parallel
  predicted <- pbmclapply(X = 1:length(models), 
                          FUN = application_predict, 
                          models = models,
                          frame = frame,
                          ncomp = ncomp,
                          mc.preschedule = TRUE, 
                          mc.set.seed = FALSE,
                          mc.cores = threads,
                          mc.cleanup = TRUE)
  
  #-----------------------------------------------------------------------------
  # Evaluate performance
  application_performance <- function(X,
                                      species_split,
                                      predicted) {
    # Predict new values
    iteration <- predicted[[X]][,-1]
    species_names <- colnames(iteration)
    predicted_species <- apply(iteration, 1, FUN = which.max)
    predicted_species <- species_names[predicted_species]
    
    # Table
    tab <- table(predicted_species, species_split)
    
    # Confusion Matrix
    cm <- confusionMatrix(tab)
    
    # Export results
    results <- data.table(iteration = X)
    results <- cbind(results, matrix(cm$overall, nrow = 1))
    colnames(results)[2:8] <- names(cm$overall)
    
    byClass <- cbind(Class = rownames(cm$byClass),
                     as.data.table(cm$byClass))
    
    results <- cbind(results, byClass)
    
    return(results)
    
  }
  
  # Performance in parallel
  performance <- pbmclapply(X = 1:length(models),
                            FUN = application_performance,
                            species_split,
                            predicted,
                            mc.preschedule = TRUE, 
                            mc.set.seed = FALSE,
                            mc.cores = threads,
                            mc.cleanup = TRUE)
  
  performance <- do.call(rbind, performance)
  performance$Class <- gsub("Class: ", "", performance$Class)
  
  #-----------------------------------------------------------------------------
  # Get mean and sd of predicted probabilities values
  
  predicted <- do.call(rbind, predicted)
  
  #-----------------------------------------------------------------------------
  # Return results
  return(list(performance = performance,
              predicted = predicted))
  
}


#-------------------------------------------------------------------------------
# Model Performance Function for LDA
model_performance_lda <- function(meta_split,
                                  species_split, 
                                  spectra_split,
                                  models,
                                  threads = 1) {
  
  # Data for predicting
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"
  
  #-----------------------------------------------------------------------------
  # Apply prediction function
  application_predict <- function(X,
                                  models,
                                  frame) {
    # Predict new values
    predicted <- predict(models[[X]], 
                         newdata = frame,
                         type = "prob")
    
    predicted <- as.data.table(cbind(iteration = X,
                                     predicted))
    
    #predicted <- as.data.table(predicted)
    
    return(predicted)
    
  }
  
  # Predicted in parallel
  predicted <- pbmclapply(X = 1:length(models), 
                          FUN = application_predict, 
                          models = models,
                          frame = frame,
                          mc.preschedule = TRUE, 
                          mc.set.seed = FALSE,
                          mc.cores = threads,
                          mc.cleanup = TRUE)
  
  #-----------------------------------------------------------------------------
  # Evaluate performance
  application_performance <- function(X,
                                      species_split,
                                      predicted) {
    # Predict new values
    iteration <- predicted[[X]][,-1]
    species_names <- colnames(iteration)
    predicted_species <- apply(iteration, 1, FUN = which.max)
    predicted_species <- species_names[predicted_species]
    
    # Table
    tab <- table(predicted_species, species_split)
    
    # Confusion Matrix
    cm <- confusionMatrix(tab)
    
    # Export results
    results <- data.table(iteration = X)
    results <- cbind(results, matrix(cm$overall, nrow = 1))
    colnames(results)[2:8] <- names(cm$overall)
    
    byClass <- cbind(Class = rownames(cm$byClass),
                     as.data.table(cm$byClass))
    
    results <- cbind(results, byClass)
    
    return(results)
    
  }
  
  # Performance in parallel
  performance <- pbmclapply(X = 1:length(models),
                            FUN = application_performance,
                            species_split,
                            predicted,
                            mc.preschedule = TRUE, 
                            mc.set.seed = FALSE,
                            mc.cores = threads,
                            mc.cleanup = TRUE)
  
  performance <- do.call(rbind, performance)
  performance$Class <- gsub("Class: ", "", performance$Class)
  
  #-----------------------------------------------------------------------------
  # Get mean and sd of predicted probabilities values
  
  predicted <- do.call(rbind, predicted)
  
  #-----------------------------------------------------------------------------
  # Return results
  return(list(performance = performance,
              predicted = predicted))
  
}

