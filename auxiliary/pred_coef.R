####################################
## TRAITS

pred_coef <- function(spectra, coef) {
  
  # Match columns between spectra and coef
  cbands <- colnames(coef)[c(-1, -2)]
  spectra_oi <- spectra[, ..cbands]
  
  if(length(cbands) != ncol(spectra_oi)) {
    stop("The number bands and coefficients needs to match")
  }
  
  # Define intercept
  coefficients <- as.matrix(coef[, c(-1, -2)])
  intercept <- coef$`(Intercept)`
  
  #Collector
  predicted_iterations <- data.table()
  
  # create progress bar
  pb <- txtProgressBar(min = 0, max = nrow(coefficients), style = 3)
  
  #Predict values
  for(ii in 1:nrow(coefficients)) {
    
    # update progress bar
    setTxtProgressBar(pb, ii)
    
    predicted <- as.matrix(spectra_oi) %*% as.numeric(coefficients[ii,])
    predicted <- predicted + intercept[ii]
    predicted_iterations <- cbind(predicted_iterations, predicted)
    
  }
  
  colnames(predicted_iterations) <- paste0("Iteration_", 1:ncol(predicted_iterations))
  
  return(predicted_iterations)

} 

####################################
### PLSDA

predict_plsda <- function(spectra, coefficients, classes, true_classes, meta, output_file = "plsda_probabilities_with_validation.csv") {
  library(data.table)
  
  # Ensure spectra columns match predictor names in coefficients
  predictors <- unique(coefficients$Predictor)
  if (!all(predictors %in% colnames(spectra))) {
    stop("Spectra data is missing some required predictor columns.")
  }
  spectra <- spectra[, ..predictors, with = FALSE]  # Match columns
  
  # Initialize a list to store predictions
  predictions <- lapply(unique(coefficients$Iteration), function(iter) {
    # Extract coefficients for the current iteration
    coef_matrix <- as.matrix(coefficients[Iteration == iter, -c("Predictor", "Iteration"), with = FALSE])
    
    # Compute PLSDA scores
    scores <- as.matrix(spectra) %*% coef_matrix
    
    # Softmax for class probabilities
    softmax <- function(x) exp(x) / rowSums(exp(x))
    probabilities <- softmax(scores)
    
    # Assign predicted classes
    predicted_classes <- classes[max.col(probabilities)]
    
    list(probabilities = probabilities, classes = predicted_classes)
  })
  
  # Aggregate results across iterations
  all_probabilities <- do.call(abind::abind, c(lapply(predictions, `[[`, "probabilities"), along = 3))  # Combine across iterations
  mean_probabilities <- apply(all_probabilities, c(1, 2), mean)  # Average probabilities
  
  # Create final predicted classes
  final_classes <- classes[max.col(mean_probabilities)]
  
  # Format probabilities for output
  probabilities_dt <- as.data.table(mean_probabilities)
  setnames(probabilities_dt, colnames(probabilities_dt), paste0("prob_", classes))  # Append "prob_" to class names
  probabilities_dt[, Sample := seq_len(nrow(probabilities_dt))]
  
  # Add true class information
  probabilities_dt[, True_Class := true_classes]
  probabilities_dt[, Predicted_Class := final_classes]
  
  # Add predicted class probabilities
  probabilities_dt[, Predicted_Class_prob := probabilities_dt[matrix(seq_len(.N), nrow = .N), paste0("prob_", Predicted_Class)]]
  
  # Add correct prediction column
  probabilities_dt[, Correct_Prediction := True_Class == Predicted_Class]
  
  # Merge with metadata
  final_output <- merge(probabilities_dt, meta, by.x = "Sample", by.y = "sample", all.x = TRUE)
  
  # Save final output to file
  fwrite(final_output, paste0(output_path, output_file))
  
  # Return validation results with metadata
  return(final_output)
}

####################################
### Generate confusion matrix based on coefficients
# Function to compute and display confusion matrix
coef_confusion_matrix <- function(predictions_output) {
  # Ensure required columns are present
  if (!all(c("True_Class", "Predicted_Class") %in% colnames(predictions_output))) {
    stop("The predictions output must contain 'True_Class' and 'Predicted_Class' columns.")
  }
  
  # Extract true and predicted classes
  true_classes <- predictions_output$True_Class
  predicted_classes <- predictions_output$Predicted_Class
  
  # Align the levels of Predicted_Class with True_Class
  true_classes <- factor(true_classes)  # Ensure true_classes is a factor
  predicted_classes <- factor(predicted_classes, levels = levels(true_classes))  # Align levels
  
  # Compute confusion matrix
  confusion <- confusionMatrix(as.factor(predicted_classes), as.factor(true_classes))
  
  # Display confusion matrix
  print(confusion$table)
  
  # Display accuracy and other metrics
  cat("\nAccuracy:", round(confusion$overall["Accuracy"], 3), "\n")
  cat("\nOther Metrics:\n")
  print(confusion$byClass)
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
