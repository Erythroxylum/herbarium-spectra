################################################################################
# Calculate optimal number of components from downsampling plsda model objects
################################################################################
getNComps = function(directory, ncomp) {
  require(matrixStats)
  models = list.files(path = directory)
  accuracies = matrix(nrow = ncomp)
  
  for(i in 1:length(models)) {
    model = readRDS(paste(directory, models[i], sep = "/"))
    accuracies = cbind(accuracies, as.matrix(model))
  }
  
  accuracies = accuracies[,-1]
  meanAcc = as.matrix(rowMeans(accuracies))
  sdAcc = as.matrix(rowSds(accuracies))
  lowerSd = meanAcc - (2 * sdAcc)
  meanAcc = meanAcc[1:which.max(meanAcc),]
  #get accuracies less than two standard deviations of the highest mean accuracy
  lowerMeans = meanAcc[meanAcc < lowerSd[which.max(meanAcc)]]
  optimalComp = length(lowerMeans) + 1 #first component with an average 
  # accuracy within 2 sd of the component with the highest mean accuracy
  
  return(optimalComp)
}

saveRDS(getNComps, "functions/getNComps.rds")
