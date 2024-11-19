# LDA VIP

lda_varImp <- function(models) {
  # Initialize an empty data frame to store variable importance for each model
  var_imp_df <- data.frame(Model = integer(), Variable = character(), Importance = numeric(), stringsAsFactors = FALSE)
  
  # Loop over each model in the list
  for (i in seq_along(models)) {
    model <- models[[i]]
    
    # Check if the model is a valid caret train object
    if (inherits(model, "train") && model$method == "lda") {
      try({
        # Calculate variable importance using varImp
        var_imp <- varImp(model, scale = TRUE)
        
        # Convert to a data frame, add model index, and format the output
        var_imp_model <- as.data.frame(var_imp$importance)
        var_imp_model$Variable <- rownames(var_imp_model)
        var_imp_model <- data.frame(Model = i, var_imp_model, row.names = NULL)
        
        # Combine with the main data frame
        var_imp_df <- rbind(var_imp_df, var_imp_model)
        
      }, silent = TRUE) # Handle any error without stopping the loop
    } else {
      warning(paste("Model", i, "is not an LDA model from caret; skipping."))
    }
  }
  
  # Check if variable importance was calculated for any model
  if (nrow(var_imp_df) == 0) {
    stop("No variable importance scores were calculated. Ensure models are LDA models from caret.")
  }
  
  return(var_imp_df)
}
