#-------------------------------------------------------------------------------
# Model Performance Function
model_performance <- function(meta_split,
                              traits_split, 
                              spectra_split,
                              models,
                              ncomp,
                              threads = 1) {
  
  # Data for predicting
  frame <- cbind(traits_split, spectra_split)
  colnames(frame)[1] <- "trait"
  
  #-----------------------------------------------------------------------------
  # Apply prediction function
  application_predict <- function(X,
                                  models,
                                  frame,
                                  ncomp) {
    # Predict new values
    predicted <- predict(models[[X]], 
                         newdata = frame,
                         ncom = ncomp)[,,1]
    
    return(predicted)
    
  }
  
  # Predicted in parallel
  predicted <- pbmclapply(X = 1:length(models), 
                          FUN = application_predict, 
                          models = models,
                          frame = frame,
                          ncomp = ncomp,
                          mc.preschedule = TRUE, 
                          mc.set.seed = FALSE,
                          mc.cores = threads,
                          mc.cleanup = TRUE)
  
  #-----------------------------------------------------------------------------
  # Evaluate performance
  application_performance <- function(X,
                                      trait_observed,
                                      predicted) {
    # Predict new values
    iteration <- predicted[[X]]
    
    # Performance at the Axes level
    per <- parameters(obs = trait_observed, 
                      pred = predicted[[X]])
    
    return(per)
    
  }
  
  # Performance in parallel
  performance <- pbmclapply(X = 1:length(models),
                            FUN = application_performance,
                            trait_observed = frame$trait,
                            predicted = predicted,
                            mc.preschedule = TRUE, 
                            mc.set.seed = FALSE,
                            mc.cores = threads,
                            mc.cleanup = TRUE)
  
  performance <- do.call(rbind, performance)
  performance <- cbind(iteration = 1:nrow(performance), 
                       nsamp = nrow(frame),
                       performance)
  
  #-----------------------------------------------------------------------------
  # Get mean and sd of predicted values and residuals
  
  # Get mean and sd residuals
  predicted <- do.call(cbind, predicted)
  observed_matrix <- matrix(rep(frame$trait, length(models)), 
                            nrow = length(frame$trait),
                            byrow = FALSE)
  residuals <- observed_matrix - predicted
  residuals_mean <- rowMeans(residuals)
  residuals_sd <- apply(residuals, 1, sd)
  
  # Get mean and sd of predicted values
  predicted_mean <- rowMeans(predicted)
  predicted_sd <- apply(predicted, 1, sd)
  
  # All together 
  predicted_summary <- meta_split[, c("accession", "leaf", "scan", "sample")]
  predicted_summary$observed <- frame$trait
  predicted_summary$predicted_mean <- predicted_mean
  predicted_summary$predicted_sd <- predicted_sd
  predicted_summary$residuals_mean <- residuals_mean
  predicted_summary$residuals_sd <- residuals_sd
  
  
  #-----------------------------------------------------------------------------
  # Return results
  return(list(performance = performance,
              predicted = predicted_summary))
  
}

#-------------------------------------------------------------------------------
# Parameters of performance to estimate 1D
parameters <- function(obs, pred) {
  
  linear <- lm(obs ~ pred, na.action = na.exclude)
  linear_summary <- summary(linear)
  
  R2 <- linear_summary$r.squared
  intercept <- linear_summary$coefficients[1,1]
  slope <- linear_summary$coefficients[2,1]
  BIAS <- mean(obs-pred, na.rm = TRUE)
  RMSE <- sqrt(sum(((obs-pred)^2))/length(obs))
  perRMSE <- sqrt(sum(((obs-pred)^2))/length(obs))/(quantile(obs, 0.99) - quantile(obs, 0.01))
  names(perRMSE) <- NULL
  
  
  return(round(data.table(R2 = R2,
                          BIAS = BIAS,
                          RMSE = RMSE,
                          perRMSE = perRMSE,
                          intercept = intercept,
                          slope = slope), 6))
  
}
