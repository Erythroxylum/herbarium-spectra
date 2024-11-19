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
  predicted_summary <- meta_split[, c("species", "genus","family","order","class","accession", "leaf", "scan", "doy", "growthForm", "absoluteAge", "herbQuality", "damage", "glue", "leafStage", "greenIndex")]
  predicted_summary$observed <- frame$trait
  predicted_summary$predicted_mean <- predicted_mean
  predicted_summary$predicted_sd <- predicted_sd
  predicted_summary$predicted_se <- predicted_se
  predicted_summary$residuals_mean <- residuals_mean
  predicted_summary$residuals_sd <- residuals_sd
  predicted_summary$residuals_se <- residuals_se
  
  # Aggregate predictions and residuals by accession
  predicted_summary_avg <- predicted_summary[, .(
    observed = mean(observed, na.rm = TRUE),
    predicted_mean = mean(predicted_mean, na.rm = TRUE),
    predicted_sd = sqrt(mean(predicted_sd^2, na.rm = TRUE)),
    predicted_rmse = sqrt(mean(predicted_se^2, na.rm = TRUE)),
    residuals_mean = mean(residuals_mean, na.rm = TRUE),
    residuals_sd = sqrt(mean(residuals_sd^2, na.rm = TRUE)),
    residuals_rmse = sqrt(mean(residuals_se^2, na.rm = TRUE)),
    species = first(species),
    genus = first(genus),
    family = first(family),
    order = first(order),
    class = first(class),
    doy = first(doy),
    growthForm = first(growthForm),
    absoluteAge = first(absoluteAge),
    herbQuality = first(herbQuality),
    damage = first(damage),
    glue = first(glue),
    leafStage = first(leafStage),
    greenIndex = first(greenIndex)
  ), by = accession]
  
  #----------------------------------------------------------------------------- 
  # Return results
  return(list(performance = performance,
              predicted = predicted_summary,
              predicted_by_accession = predicted_summary_avg))
}







#-----------------------------------------------------------------------------
# Kothari dataset
#-----------------------------------------------------------------------------
model_performance_kothari <- function(meta_split,
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
  
  # Residuals statistics
  residuals_mean <- rowMeans(residuals)
  residuals_sd <- apply(residuals, 1, sd)
  residuals_se <- residuals_sd / sqrt(length(models))
  
  # Predictions statistics
  predicted_mean <- rowMeans(predicted)
  predicted_sd <- apply(predicted, 1, sd)
  predicted_se <- predicted_sd / sqrt(length(models))
  
  # Combine with meta information and summarize by `accession`
  predicted_summary <- meta_split[, c("species", "genus","accession","growthForm","Project", 
                                      "greenIndex","Discoloration")]
  predicted_summary$observed <- frame$trait
  predicted_summary$predicted_mean <- predicted_mean
  predicted_summary$predicted_sd <- predicted_sd
  predicted_summary$predicted_se <- predicted_se
  predicted_summary$residuals_mean <- residuals_mean
  predicted_summary$residuals_sd <- residuals_sd
  predicted_summary$residuals_se <- residuals_se
  
  # Aggregate predictions and residuals by accession
  predicted_summary_avg <- predicted_summary[, .(observed = mean(observed, na.rm = TRUE),
                                                 predicted_mean = mean(predicted_mean, na.rm = TRUE),
                                                 predicted_sd = sqrt(mean(predicted_sd^2, na.rm = TRUE)),
                                                 predicted_rmse = sqrt(mean(predicted_se^2, na.rm = TRUE)),
                                                 residuals_mean = mean(residuals_mean, na.rm = TRUE),
                                                 residuals_sd = sqrt(mean(residuals_sd^2, na.rm = TRUE)),
                                                 residuals_rmse = sqrt(mean(residuals_se^2, na.rm = TRUE)),
                                                 species = first(species),
                                                 genus = first(genus),
                                                 growthForm = first(growthForm),
                                                 Project = first(Project),
                                                 Discoloration = first(Discoloration),
                                                 greenIndex = first(greenIndex)), 
                                             by = accession]
  
  #----------------------------------------------------------------------------- 
  # Return results
  return(list(performance = performance,
              predicted = predicted_summary,
              predicted_by_accession = predicted_summary_avg))
  
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
  
  
  return(round(data.table(R2 = R2,
                          BIAS = BIAS,
                          RMSE = RMSE,
                          perRMSE = perRMSE,
                          intercept = intercept,
                          slope = slope), 6))
  
}
