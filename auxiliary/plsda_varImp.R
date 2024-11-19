# Define the function to extract variable importance from all PLSDA models
plsda_varImp <- function(models) {
  # Check if the input is a list of models
  if (!is.list(models)) {
    stop("Error: Input is not a list of models.")
  }
  
  # Initialize an empty data frame to store variable importance scores for each model
  var_imp_df <- data.frame(Model = integer(), Variable = character(), Importance = numeric(), stringsAsFactors = FALSE)
  
  # Iterate through each model in the list
  for (i in seq_along(models)) {
    model <- models[[i]]
    
    # Check that the model contains a 'finalModel' and apply varImp
    if ("finalModel" %in% names(model)) {
      try({
        # Calculate variable importance
        var_imp <- varImp(model, scale = TRUE)
        
        # Convert to a data frame and add model information
        var_imp_model <- as.data.frame(var_imp$importance)
        var_imp_model$Variable <- rownames(var_imp_model)
        var_imp_model <- data.frame(Model = i, var_imp_model, row.names = NULL)
        
        # Store the importance data
        var_imp_df <- rbind(var_imp_df, var_imp_model)
        
      }, silent = TRUE) # Suppress errors and continue with next model
    } else {
      warning(paste("Model", i, "does not have a 'finalModel' component; skipping."))
    }
  }
  
  # Check if variable importance scores were successfully calculated
  if (nrow(var_imp_df) == 0) {
    stop("No variable importance scores were calculated. Check model compatibility with varImp().")
  }
  
  return(var_imp_df)
}





