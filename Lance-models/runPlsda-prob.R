################################################################################
# PLS function
################################################################################

runPlsda = function(iteration, spectra, className, ncomp, resampling, include_age,
                    modelDirectory, accuracyDirectory, cmDirectory,
                    probDirectory, saveModelObject) {
  require(spectrolab)
  require(caret)
  require(dplyr)
  require(rlist)
  require(matrixStats)
  require(mlbench)
  require(MASS)
  require(parallel)
  require(CAST)

  # Load and preprocess data
  spec_raw = spectra
  spec_mat = as.matrix(spectra)
  spec_df1 = as.data.frame(spec_raw)
  spec_df2 = as.data.frame(spec_mat)
  spec_df2 = cbind(spec_df2, spec_df1[className])
  colnames(spec_df2)[colnames(spec_df2) == className] <- className

  if (include_age == TRUE) {
    spec_df2$age = spec_df1$age
  }

  # Create partitions for training/testing using CAST
  folds <- CreateSpacetimeFolds(data = spec_df2, spacevar = "accession", k = 1)
  inTrain <- folds$index[[1]]
  inTest <- folds$indexOut[[1]]

  training <- spec_df2[inTrain, ]
  testing  <- spec_df2[inTest, ]

  # Set up model tuning
  ctrl <- trainControl(method = "repeatedcv", number = 10, sampling = resampling, repeats = 1)

  # Fit the PLS model
  plsFit <- train(as.formula(paste(className, "~ .")),
                  data = training,
                  method = "pls",
                  trControl = ctrl,
                  tuneLength = ncomp)

  # Save the PLS model
  fileName = paste(paste('pls', className, "model", age, resampling, toString(iteration), sep = "_"), "rds", sep = ".")
  saveRDS(plsFit, paste(modelDirectory, fileName, sep = "/"))

  # Predict the class probabilities on the test set
  plsProbs = predict(plsFit, newdata = testing, type = "prob")
  
  # Print the probabilities (optional)
  print(paste("Class probabilities for iteration", iteration, ":"))
  #print(plsProbs)

  # Save the probabilities to a CSV file
  probFileName = paste(paste('probabilities', className, age, resampling, toString(iteration), sep = "_"), "csv", sep = ".")
  write.csv(plsProbs, file = paste(probDirectory, probFileName, sep = "/"), row.names = FALSE)

  if (resampling == 'up') {
    # Predict the final class (based on highest probability)
    plsClasses = predict(plsFit, newdata = testing)

    # Generate the confusion matrix
    cm = confusionMatrix(data = plsClasses, as.factor(testing[[className]]))

    # Save the confusion matrix
    cmFileName = paste(paste('cm', toString(iteration), sep="_"), "rds", sep = ".")
    saveRDS(cm, paste(cmDirectory, cmFileName, sep = "/"))
  }
}
saveRDS(runPlsda, "functions/runPlsda-prob.rds")












