pred_coef <- function(spectra, coef, bands) {
  
  # Ensure bands are characters (for matching column names)
  bands <- as.character(bands)
  
  # Filter coef to only include columns matching bands
  coef_filtered <- coef[, c("(Intercept)", bands), with = FALSE]
  
  # Match columns between spectra and coef
  cbands <- intersect(colnames(spectra), bands) # Columns to use from spectra
  spectra_oi <- spectra[, ..cbands, with = FALSE] # Subset spectra
  
  # Check for mismatch between bands in spectra and coef
  if (length(cbands) != length(bands)) {
    warning("Some bands were not found in spectra or coef; filtering to match.")
  }
  
  if (length(cbands) != ncol(coef_filtered) - 1) { # Subtract 1 for intercept
    stop(paste0("Mismatch: Found ", length(cbands), 
                " bands in spectra but expected ", ncol(coef_filtered) - 1, 
                " coefficients after filtering."))
  }
  
  # Define intercept and coefficients
  coefficients <- as.matrix(coef_filtered[, -1, with = FALSE]) # Exclude intercept
  intercept <- coef_filtered$`(Intercept)`
  
  # Collector for predictions
  predicted_iterations <- data.table()
  
  # Create progress bar
  pb <- txtProgressBar(min = 0, max = nrow(coefficients), style = 3)
  
  # Predict values
  for (ii in 1:nrow(coefficients)) {
    
    # Update progress bar
    setTxtProgressBar(pb, ii)
    
    # Perform prediction
    predicted <- as.matrix(spectra_oi) %*% as.numeric(coefficients[ii, ])
    predicted <- predicted + intercept[ii]
    predicted_iterations <- cbind(predicted_iterations, predicted)
  }
  
  close(pb) # Close progress bar
  
  # Assign column names
  colnames(predicted_iterations) <- paste0("Iteration_", 1:ncol(predicted_iterations))
  
  return(predicted_iterations)
}
