#-------------------------------------------------------------------------------
# Function to get the confusion matrices

confusionMatrices_plsda <- function(meta_split,
                              species_split, 
                              spectra_split,
                              models = models_plsda,
                              ncomp = ncomp,
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
    cm <- cm$table
    
    return(cm)
    
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
  
  rownames_cm <- rownames(performance[[1]])
  colnames_cm <- colnames(performance[[1]])
  
  performance <- array(unlist(performance), dim = c(dim(performance[[1]]),length(performance)))
  
  #-----------------------------------------------------------------------------
  # Get mean and sd of the confusion values
  
  dims_values <- dim(performance)
  collector_mean <- matrix(NA, nrow = dims_values[1], ncol = dims_values[2])
  collector_sd <- matrix(NA, nrow = dims_values[1], ncol = dims_values[2])
  
  for(i in 1:dims_values[1]) {
    for(ii in 1:dims_values[2]) {
      
      collector_mean[i, ii] <- mean(performance[i,ii,])
      collector_sd[i, ii] <- sd(performance[i,ii,]) 
      
    }
  }
  
  colnames(collector_mean) <-  colnames_cm
  colnames(collector_sd) <-  colnames_cm
  rownames(collector_mean) <-  rownames_cm
  rownames(collector_sd) <-  rownames_cm
  
  
  #-----------------------------------------------------------------------------
  # Return results
  return(list(confusionMatrices_mean = collector_mean,
              confusionMatrices_sd = collector_sd))
  
  
}


confusionMatrices_lda <- function(meta_split,
                                    species_split, 
                                    spectra_split,
                                    models = models_plsda,
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
    cm <- cm$table
    
    return(cm)
    
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
  
  rownames_cm <- rownames(performance[[1]])
  colnames_cm <- colnames(performance[[1]])
  
  performance <- array(unlist(performance), dim = c(dim(performance[[1]]),length(performance)))
  
  #-----------------------------------------------------------------------------
  # Get mean and sd of the confusion values
  
  dims_values <- dim(performance)
  collector_mean <- matrix(NA, nrow = dims_values[1], ncol = dims_values[2])
  collector_sd <- matrix(NA, nrow = dims_values[1], ncol = dims_values[2])
  
  for(i in 1:dims_values[1]) {
    for(ii in 1:dims_values[2]) {
      
      collector_mean[i, ii] <- mean(performance[i,ii,])
      collector_sd[i, ii] <- sd(performance[i,ii,]) 
      
    }
  }
  
  colnames(collector_mean) <-  colnames_cm
  colnames(collector_sd) <-  colnames_cm
  rownames(collector_mean) <-  rownames_cm
  rownames(collector_sd) <-  rownames_cm
  
  
  #-----------------------------------------------------------------------------
  # Return results
  return(list(confusionMatrices_mean = collector_mean,
              confusionMatrices_sd = collector_sd))
  
  
}
