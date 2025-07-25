####################################
### PLSDA

predict_plsda <- function(coefficients_list, spectra, meta, traits) {
  # Ensure coefficients, xmeans, and ymeans are data.tables
  coefficients_dt <- as.data.table(coefficients_list$coefficients)
  xmeans_dt <- as.data.table(coefficients_list$xmeans)
  ymeans_dt <- as.data.table(coefficients_list$ymeans)
  
  # Softmax function
  softmax <- function(x) exp(x) / rowSums(exp(x))
  
  # Convert spectra to numeric matrix
  spectra <- as.matrix(spectra)
  
  # Predict probabilities across iterations
  probabilities_list <- lapply(unique(coefficients_dt$Iteration), function(iter) {
    # Extract coefficients and means for iteration
    coeff_iter <- coefficients_dt[Iteration == iter, !c("Iteration", "Predictor"), with = FALSE]
    coeff_matrix <- as.matrix(coeff_iter)
    xmeans_iter <- as.numeric(xmeans_dt[Iteration == iter, Xmean])
    ymeans_iter <- as.numeric(ymeans_dt[Iteration == iter, Ymean])
    
    # Compute intercept
    int0 <- ymeans_iter - xmeans_iter %*% coeff_matrix
    
    # Raw predictions
    raw_predictions_iter <- spectra %*% coeff_matrix + rep(int0, each = nrow(spectra))
    
    # Convert to probabilities
    probabilities_iter <- softmax(raw_predictions_iter)
    
    # Output as data.table
    probabilities_dt <- as.data.table(probabilities_iter)
    colnames(probabilities_dt) <- setdiff(colnames(coeff_iter), "Predictor")
    probabilities_dt[, Iteration := iter]
    probabilities_dt
  })
  
  # Combine across iterations
  combined_probabilities <- rbindlist(probabilities_list)
  combined_probabilities_no_iter <- combined_probabilities[, !"Iteration", with = FALSE]
  
  # Reshape into 3D array: obs × iterations × classes
  n_observations <- nrow(combined_probabilities) / length(unique(combined_probabilities$Iteration))
  n_iterations <- length(unique(combined_probabilities$Iteration))
  n_classes <- ncol(combined_probabilities_no_iter)
  
  probabilities_array <- array(as.matrix(combined_probabilities_no_iter),
                               dim = c(n_observations, n_iterations, n_classes))
  
  # Average over iterations
  average_probabilities <- apply(probabilities_array, c(1, 3), mean)
  average_probabilities_dt <- as.data.table(average_probabilities)
  setnames(average_probabilities_dt, colnames(combined_probabilities_no_iter))
  
  # Predicted class = max prob
  class_labels <- colnames(average_probabilities_dt)
  predicted_classes <- factor(class_labels[max.col(average_probabilities, ties.method = "first")],
                              levels = class_labels)
  
  # True classes (underscored)
  true_classes <- factor(sub(" ", "_", meta$scientificName), levels = class_labels)
  
  # Align factor levels before comparison
  all_levels <- union(levels(predicted_classes), levels(true_classes))
  predicted_classes <- factor(predicted_classes, levels = all_levels)
  true_classes <- factor(true_classes, levels = all_levels)
  
  # Logical vector for prediction correctness
  prediction_correct <- predicted_classes == true_classes
  
  # Add predicted class column for prob lookup
  average_probabilities_dt[, predicted_class := predicted_classes]
  
  # Add prob_predicted column using mapply
  average_probabilities_dt[, prob_predicted := mapply(function(class, row)
    row[[as.character(class)]],
    predicted_class, split(.SD, seq_len(.N))),
    .SDcols = class_labels]
  
  # Final output
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
