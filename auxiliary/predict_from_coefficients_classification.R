####################################
### PLSDA

predict_plsda <- function(coefficients_list, spectra, meta, traits) {
  # Ensure coefficients, xmeans, and ymeans are treated as data.tables
  coefficients_dt <- as.data.table(coefficients_list$coefficients)
  xmeans_dt <- as.data.table(coefficients_list$xmeans)
  ymeans_dt <- as.data.table(coefficients_list$ymeans)
  
  # Define softmax function
  softmax <- function(x) exp(x) / rowSums(exp(x))
  
  #make spectra into matrix numeric
  spectra <- as.matrix(spectra)
  
  # Calculate class probabilities across iterations
  probabilities_list <- lapply(unique(coefficients_dt$Iteration), function(iter) {
    # Extract coefficients, xmeans, and ymeans for the current iteration
    coeff_iter <- coefficients_dt[Iteration == iter, !c("Iteration", "Predictor"), with = FALSE]
    coeff_matrix <- as.matrix(coeff_iter)
    xmeans_iter <- as.numeric(xmeans_dt[Iteration == iter, Xmean])
    ymeans_iter <- as.numeric(ymeans_dt[Iteration == iter, Ymean])
    
    # Compute the intercept (int0)
    int0 <- ymeans_iter - xmeans_iter %*% coeff_matrix
    
    # Compute raw predictions
    raw_predictions_iter <- spectra %*% coeff_matrix + rep(int0, each = nrow(spectra))
    
    # Compute probabilities
    probabilities_iter <- softmax(raw_predictions_iter)
    
    # Return probabilities as a data.table
    probabilities_dt <- as.data.table(probabilities_iter)
    colnames(probabilities_dt) <- setdiff(colnames(coeff_iter), "Predictor")
    probabilities_dt[, Iteration := iter]
    probabilities_dt
  })
  
  # Combine probabilities across all iterations
  combined_probabilities <- rbindlist(probabilities_list)
  
  # Remove the "Iteration" column
  combined_probabilities_no_iter <- combined_probabilities[, !"Iteration", with = FALSE]
  
  # Reshape to 3D array for easier averaging
  n_observations <- nrow(combined_probabilities) / length(unique(combined_probabilities$Iteration))
  n_iterations <- length(unique(combined_probabilities$Iteration))
  probabilities_array <- array(as.matrix(combined_probabilities_no_iter), 
                               dim = c(n_observations, n_iterations, ncol(combined_probabilities_no_iter)))
  
  # Average across iterations
  average_probabilities <- apply(probabilities_array, c(1, 3), mean)
  
  # Convert back to a data.table
  average_probabilities_dt <- as.data.table(average_probabilities)
  setnames(average_probabilities_dt, colnames(combined_probabilities_no_iter))  # Match column names
  
  # Generate predicted classes
  classes <- colnames(combined_probabilities_no_iter)  # Class names
  predicted_classes <- classes[max.col(average_probabilities)]  # Get class with highest probability
  predicted_classes <- as.factor(predicted_classes)
  
  # True classes (convert to factors)
  true_classes <- as.factor(meta$scientificName)
  true_classes <- as.factor(sub(" ", "_", true_classes))
  
  # Create a column for whether the prediction is correct
  prediction_correct <- predicted_classes == true_classes
  
  # generate column for the probability of the predicted class
  average_probabilities_dt[, prob_predicted := mapply(function(class, row) row[[as.character(class)]],
                                               predicted_class, split(.SD, seq_len(.N)), 
                                               SIMPLIFY = TRUE), 
                    .SDcols = colnames(predictions_table)[sapply(colnames(predictions_table), function(x) is.numeric(predictions_table[[x]]))]]
  
  # Assemble the final output data.table
  output_dt <- data.table(
    true_class = true_classes,
    predicted_class = predicted_classes,
    correct = prediction_correct,
    average_probabilities_dt,
    meta,
    traits
  )
  
  # Return the final output
  return(output_dt)
}


predict_plsda <- function(coefficients_list, spectra, meta, traits) {
  # Ensure coefficients, xmeans, and ymeans are treated as data.tables
  coefficients_dt <- as.data.table(coefficients_list$coefficients)
  xmeans_dt <- as.data.table(coefficients_list$xmeans)
  ymeans_dt <- as.data.table(coefficients_list$ymeans)
  
  # Define softmax function
  softmax <- function(x) exp(x) / rowSums(exp(x))
  
  # Make spectra into matrix numeric
  spectra <- as.matrix(spectra)
  
  # Calculate class probabilities across iterations
  probabilities_list <- lapply(unique(coefficients_dt$Iteration), function(iter) {
    # Extract coefficients, xmeans, and ymeans for the current iteration
    coeff_iter <- coefficients_dt[Iteration == iter, !c("Iteration", "Predictor"), with = FALSE]
    coeff_matrix <- as.matrix(coeff_iter)
    xmeans_iter <- as.numeric(xmeans_dt[Iteration == iter, Xmean])
    ymeans_iter <- as.numeric(ymeans_dt[Iteration == iter, Ymean])
    
    # Compute the intercept (int0)
    int0 <- ymeans_iter - xmeans_iter %*% coeff_matrix
    
    # Compute raw predictions
    raw_predictions_iter <- spectra %*% coeff_matrix + rep(int0, each = nrow(spectra))
    
    # Compute probabilities
    probabilities_iter <- softmax(raw_predictions_iter)
    
    # Return probabilities as a data.table
    probabilities_dt <- as.data.table(probabilities_iter)
    colnames(probabilities_dt) <- setdiff(colnames(coeff_iter), "Predictor")
    probabilities_dt[, Iteration := iter]
    probabilities_dt
  })
  
  # Combine probabilities across all iterations
  combined_probabilities <- rbindlist(probabilities_list)
  
  # Remove the "Iteration" column
  combined_probabilities_no_iter <- combined_probabilities[, !"Iteration", with = FALSE]
  
  # Reshape to 3D array for easier averaging
  n_observations <- nrow(combined_probabilities) / length(unique(combined_probabilities$Iteration))
  n_iterations <- length(unique(combined_probabilities$Iteration))
  probabilities_array <- array(as.matrix(combined_probabilities_no_iter), 
                               dim = c(n_observations, n_iterations, ncol(combined_probabilities_no_iter)))
  
  # Average across iterations
  average_probabilities <- apply(probabilities_array, c(1, 3), mean)
  
  # Convert back to a data.table
  average_probabilities_dt <- as.data.table(average_probabilities)
  setnames(average_probabilities_dt, colnames(combined_probabilities_no_iter))  # Match column names
  
  # Generate predicted classes
  classes <- colnames(combined_probabilities_no_iter)  # Class names
  predicted_classes <- classes[max.col(average_probabilities)]  # Get class with highest probability
  predicted_classes <- as.factor(predicted_classes)
  
  # True classes (convert to factors)
  true_classes <- as.factor(meta$species)
  true_classes <- as.factor(sub(" ", "_", true_classes))
  
  # Create a column for whether the prediction is correct
  prediction_correct <- predicted_classes == true_classes
  
  # Generate column for the probability of the predicted class
  average_probabilities_dt[, prob_predicted := mapply(function(class, row) row[[as.character(class)]],
                                                      predicted_classes, split(.SD, seq_len(.N)), 
                                                      SIMPLIFY = TRUE), 
                           .SDcols = colnames(average_probabilities_dt)]
  
  # Assemble the final output data.table
  output_dt <- data.table(
    true_class = true_classes,
    predicted_class = predicted_classes,
    correct = prediction_correct,
    average_probabilities_dt,
    meta,
    traits
  )
  
  return(output_dt)
}



####################################
### LDA # not functioning. Not sure how LDs or scaling values translate to coefficients

predict_lda <- function(spectra, coefficients, num_lds, classes) {
  # Ensure spectra columns match predictor names in coefficients
  predictors <- unique(coefficients$Predictor)
  spectra <- spectra[, ..predictors, with = FALSE]  # Match columns
  
  # Extract LD coefficients
  ld_coefficients <- as.matrix(coefficients[, 1:num_lds, with = FALSE])
  
  # Compute LD values
  ld_values <- as.matrix(spectra) %*% t(ld_coefficients)
  
  # Assign classes based on the LD values
  # The order of classes corresponds to the order of LDs
  predicted_classes <- apply(ld_values, 1, function(x) classes[which.max(x)])
  
  return(predicted_classes)
}
