# Define a function to calculate variable importance for a list of models using caret's varImp()
# Define the function to extract variable importance from all LDA models
lda_varImp <- function(models_lda) {
  # Initialize an empty list to store variable importance data for all models
  lda_varImp_all_models <- list()
  
  # Loop over each model in the models_lda list
  for (i in 1:length(models_lda)) {
    # Access the final model of the current LDA model
    lda_model <- models_lda[[i]]$finalModel
    
    # Extract the variable importance (scaling matrix)
    lda_varImp <- lda_model$scaling
    
    # Convert the scaling matrix to a data frame for easier interpretation
    lda_varImp_df <- as.data.frame(lda_varImp)
    
    # Set the row and column names to match your feature names (if necessary)
    colnames(lda_varImp_df) <- paste0("LD", 1:ncol(lda_varImp_df))
    
    # Add a column for the variable names
    lda_varImp_df$Variable <- rownames(lda_varImp_df)
    
    # Store the variable importance data for the current model
    lda_varImp_all_models[[i]] <- lda_varImp_df
  }
  
  # Combine all variable importance data into a single data frame
  lda_varImp_combined <- do.call(rbind, lda_varImp_all_models)
  
  # Optionally, add a column indicating which model the variable importance is from
  lda_varImp_combined$Model <- rep(1:length(models_lda), each = nrow(lda_varImp_combined) / length(models_lda))
  
  # Return the combined variable importance data
  return(lda_varImp_combined)
}






