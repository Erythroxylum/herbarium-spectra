################################################################################
# Calculate optimal number of components from downsampling plsda model objects
################################################################################

# Revised getNComps function using selectNcomp with one-minus-sigma rule
getNComps = function(directory, ncomp) {
  require(matrixStats)
  require(pls)  # ensure pls package is loaded for selectNcomp()
  
  models = list.files(path = directory)
  accuracies = matrix(nrow = ncomp)
  
  for (i in 1:length(models)) {
    model = readRDS(paste(directory, models[i], sep = "/"))
    accuracies = cbind(accuracies, as.matrix(model))
  }
  
  accuracies = accuracies[, -1]  # Remove the first empty column
  
  # Calculate the mean and standard deviation of accuracies across all models
  meanAcc = as.matrix(rowMeans(accuracies))
  sdAcc = as.matrix(rowSds(accuracies))
  
  # Fit a partial least squares model to obtain CV error for components
  plsFit <- plsr(as.matrix(accuracies) ~ 1:ncomp, validation = "CV")
  
  # Select the optimal number of components using the one-minus-sigma rule
  optimalComp = selectNcomp(plsFit, method = "onesigma", ncomp = ncomp)
  
  return(optimalComp)
}

saveRDS(getNComps, "functions/getNComps_cv.rds")