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
