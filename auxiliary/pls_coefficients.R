#-------------------------------------------------------------------------------
# Extract the pls coefficients

pls_coefficients <- function(models = models,
                             ncomp) {
  
  # Collector
  coefficients <- data.table()
  
  # Number of models
  iterations <- length(models)
  
  # Coefficients
  for(i in 1:iterations) {
    
    # Select model
    object <- models[[i]]
    
    # Get the coefficients
    coef <- coef(object, ncomp = ncomp, intercept=TRUE)[,,1]
    coefficients <- rbind(coefficients, matrix(coef, nrow = 1))
    
  }
  
  colnames(coefficients) <- names(coef(models[[1]], ncomp = ncomp, intercept=TRUE)[,,1])
  coefficients <- cbind(iteration = 1:nrow(coefficients),
                        coefficients)
  
  
  return(coefficients)
  
}

#-------------------------------------------------------------------------------
# Extract the plsda coefficients

plsda_coefficients <- function(models, ncomp) {
  # Number of models (iterations)
  iterations <- length(models)
  
  # Extract coefficients for each model
  coefficients <- lapply(1:iterations, function(i) {
    # Extract coefficients for the desired number of components
    coef_matrix <- models[[i]]$finalModel$coefficients[, , ncomp]
    
    # Convert to a data.table
    coef_dt <- as.data.table(as.matrix(coef_matrix))
    coef_dt[, Predictor := rownames(coef_matrix)]  # Add predictor names as a column
    coef_dt[, Iteration := i]  # Add iteration identifier
    
    return(coef_dt)
  })
  
  # Combine all iterations into one data.table
  coefficients_combined <- rbindlist(coefficients, use.names = TRUE, fill = TRUE)
  
  # Return the combined coefficients table
  return(coefficients_combined)
}

#-------------------------------------------------------------------------------
# Extract the LDA coefficients

lda_coefficients <- function(models) {
  # Number of models (iterations)
  iterations <- length(models)
  
  # Extract coefficients for each model
  coefficients <- lapply(1:iterations, function(i) {
    # Extract scaling coefficients for the current model
    coef_matrix <- models[[i]]$finalModel$scaling  # Extract LDA scaling coefficients
    
    # Convert scaling coefficients to a data.table
    coef_dt <- as.data.table(as.matrix(coef_matrix))  # Convert to a data.table
    coef_dt[, Predictor := rownames(coef_matrix)]  # Add predictor names as a column
    coef_dt[, Iteration := i]  # Add iteration identifier
    
    return(coef_dt)
  })
  
  # Combine all iterations into one data.table
  coefficients_combined <- rbindlist(coefficients, use.names = TRUE, fill = TRUE)
  
  # Return the combined coefficients table
  return(coefficients_combined)
}
