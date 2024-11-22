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
  
  # Performance in parallel
  performance <- pbmclapply(X = 1:length(models),
                            FUN = application_performance,
                            meta_split,
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
  # Get mean, sd, and standard error of predicted values and residuals
  
  # Combine predictions from iterations
  predicted <- do.call(cbind, predicted)
  observed_matrix <- matrix(rep(frame$trait, length(models)), 
                            nrow = length(frame$trait),
                            byrow = FALSE)
  residuals <- observed_matrix - predicted
  
  # Residuals statistics
  residuals_mean <- rowMeans(residuals)
  residuals_sd <- apply(residuals, 1, sd)
  residuals_se <- residuals_sd / sqrt(length(models))
  
  # Predictions statistics
  predicted_mean <- rowMeans(predicted)
  predicted_sd <- apply(predicted, 1, sd)
  predicted_se <- predicted_sd / sqrt(length(models))
  
  # Combine with meta information and summarize by `accession`
  predicted_summary <- meta_split
  predicted_summary$observed <- frame$trait
  predicted_summary$predicted_mean <- predicted_mean
  predicted_summary$predicted_sd <- predicted_sd
  predicted_summary$predicted_se <- predicted_se
  predicted_summary$residuals_mean <- residuals_mean
  predicted_summary$residuals_sd <- residuals_sd
  predicted_summary$residuals_se <- residuals_se
  
  #----------------------------------------------------------------------------- 
  # Return results
  return(list(performance = performance,
              predicted = predicted_summary))
}

# ------------------------------------------------------------------------------
# Application performance
application_performance <- function(X,
                                    meta_split,
                                    trait_observed,
                                    predicted) {
  
  # Ave
  ave <- data.table(accession = meta_split$accession,
                    trait_observed = trait_observed,
                    predicted = predicted[[X]])
  ave <- ave[, .(mean(trait_observed), mean(predicted)), 
             by = "accession"]
  colnames(ave) <- c("accession", "trait_observed", "predicted")
  
  # Performance at the Axes level
  per <- parameters(obs = ave$trait_observed, 
                    pred = ave$predicted)
  
  return(per)
  
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
  perRMSE <- sqrt(sum(((obs-pred)^2))/length(obs))/(quantile(obs, 0.99, na.rm = TRUE) - quantile(obs, 0.01, na.rm = TRUE))
  names(perRMSE) <- NULL
  
  
  return(round(data.table(n = length(obs),
                          R2 = R2,
                          BIAS = BIAS,
                          RMSE = RMSE,
                          perRMSE = perRMSE,
                          intercept = intercept,
                          slope = slope), 6))
  
}
