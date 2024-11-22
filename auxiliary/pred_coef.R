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
