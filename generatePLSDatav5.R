################################################################################
# Revised plsda workflow
################################################################################

generatePLSData = function(spectra, className, includeAge, ncomps,
                           numIterations, baseDirectory) {

  #############
  # Setup
  #############
  
  require(parallel)
  
  # functions
  runPlsda_path = "functions/runPlsda-prob.rds"
  getNComps_path = "functions/getNComps_cv.rds"
  getData_path = "functions/getData.rds"
  plotConfusionMatrix_path = "functions/plotConfusionMatrix.rds"
  plotVip_path = "functions/plotVip.rds"

  print(paste("Loading function from:", runPlsda_path))
  runPlsda = readRDS(runPlsda_path)
  
  print(paste("Loading function from:", getNComps_path))
  getNComps = readRDS(getNComps_path)
  
  print(paste("Loading function from:", getData_path))
  getData = readRDS(getData_path)
  
  print(paste("Loading function from:", plotConfusionMatrix_path))
  plotConfusionMatrix = readRDS(plotConfusionMatrix_path)
  
  print(paste("Loading function from:", plotVip_path))
  plotVip = readRDS(plotVip_path)
  
  # variables 
  if (includeAge) {
    age = 'with-age'
  } else {
    age = 'no-age'
  }
  upSamplingDirectory = paste(baseDirectory, className, age, 'upsampled-models', sep = '/' )
  accuracyDirectory = paste(baseDirectory, className, age, 'accuracies', sep = '/')
  cmDirectory = paste(baseDirectory, className, age, 'confusion-matrices', sep = '/')
  metricsDirectory = paste(baseDirectory, className, age, 'metrics', sep = '/')
  cmFinalDirectory = paste(baseDirectory, className, age, 'cm-final', sep = '/')
  vipDirectory = paste(baseDirectory, className, age, 'variable-importance', sep = '/')
  
  # Print directories for debugging
  print(paste("Upsampling Directory:", upSamplingDirectory))
  print(paste("Accuracy Directory:", accuracyDirectory))
  print(paste("Confusion Matrix Directory:", cmDirectory))
  print(paste("Metrics Directory:", metricsDirectory))
  print(paste("Final Confusion Matrix Directory:", cmFinalDirectory))
  print(paste("Variable Importance Directory:", vipDirectory))
  
  #################################
  # Classification - parallelized
  #################################
  
  # choose number of cores - this code chooses half the cores available on your
  # machine
  numCores = floor(detectCores())
  cluster = makeCluster(numCores)
  iterations = seq(1, numIterations)
  
  # this function will run PLSDA and save model objects into the specified directory
  print("Starting PLSDA with downsampling...")
  parLapply(cl = cluster, iterations, runPlsda, spectra = spectra,
            className = className, ncomp = ncomps, resampling = "down",
            include_age = includeAge, modelDirectory = '', saveModelObject = FALSE,
            cmDirectory = "", accuracyDirectory = accuracyDirectory)
  
  # get optimal number of components from downsampled model objects
  print("Determining optimal number of components...")
  optimalNumComps = getNComps(accuracyDirectory, ncomps)
  
  # run plsda with upsampling with optimal number of components
  print(paste("Starting PLSDA with upsampling using", optimalNumComps, "components..."))
  parLapply(cl = cluster, iterations, runPlsda, spectra = spectra,
            className = className, ncomp = optimalNumComps, resampling = "up",
            include_age = includeAge, modelDirectory = upSamplingDirectory, 
            cmDirectory = cmDirectory, saveModelObject = TRUE)
  
  print("Completed PLSDA.")

  # Now read and save the probabilities
  for (iteration in iterations) {
    probFileName = paste(paste('probabilities', className, age, 'up', toString(iteration), sep = "_"), "csv", sep = ".")
    probPath = paste(upSamplingDirectory, probFileName, sep = "/")
    
    if (file.exists(probPath)) {
      probs = readRDS(probPath)
      # Save the output (for example, as an RDS file)
      saveRDS(probs, file = paste(upSamplingDirectory, paste('probs', className, age, 'up', toString(iteration), "rds", sep="_"), sep="/"))
    } else {
      print(paste("Probabilities file not found for iteration", iteration))
    }
  }

  ##############################################################################
  # Get data and plot confusion matrices and variable importance values
  ##############################################################################
  
  # get overall accuracy, mean confusion matrix, and standard deviation (2) 
  # confusion matrix
  print("Collecting data for metrics and plots...")
  data = getData(directory = cmDirectory, metricsDirectory = metricsDirectory,
                 className = className, includesAge = includeAge)
  
  # Set the correct number of components in data
  data$ncomps <- optimalNumComps
  
  # Now use optimalNumComps for the matrix name
  matrixName = paste(paste(className, age, optimalNumComps, 'comps', sep = '_'),
                     'jpeg', sep = '.')
  
  # plot confusion matrix as high resolution jpeg
  plotConfusionMatrix(data$cmMean, directory = cmFinalDirectory,
                      fileName = matrixName)
  
  # plot top 5 and bottom 5 variable importance values
  plotVip(modelDirectory = upSamplingDirectory, saveDirectory = vipDirectory,
          baseFileName = paste(className, 'vip', sep = '_'))
  
  print("All tasks completed successfully!")
}

saveRDS(generatePLSData, 'functions/generatePLSDatav5.rds')
