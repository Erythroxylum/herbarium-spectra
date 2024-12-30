## PLS Variable Importance in Projection

pls_vip <- function(models = models,
                    ncomp = ncomp) {
  
  #Collector
  vip <- data.table()
  
  #Number of models
  iterations <- length(models)
  
  # VIP
  for(i in 1:iterations) {
    
    #Select model
    object <- models[[i]]
    
    #Get VIP
    vip_model <- VIP(object, ncomp)
    vip <- rbind(vip, matrix(vip_model, nrow = 1))
    
  }
  
  colnames(vip) <- names(VIP(models[[1]], ncomp))
  vip <- cbind(iteration = 1:nrow(vip),
               vip)
  
  return(vip)
  
}


##############################
## PLSDA Variable Importance in Projection

# Define the function to extract variable importance from all PLSDA models
plsda_vip <- function(models) {
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


#########################
## LDA VIP

lda_vip <- function(models_lda) {
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
