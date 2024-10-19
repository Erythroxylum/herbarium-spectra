################################################################################
# Revised plsda workflow
################################################################################

generatePLSData = function(spectra, className, includeAge, ncomps,
                           numIterations, baseDirectory) {

  #############
  # Setup
  #############
  
  require(parallel)
  
  #functions
  runPlsda = readRDS("functions/runPlsda.rds")
  getNComps = readRDS("functions/getNComps.rds")
  getData = readRDS("functions/getData.rds")
  plotConfusionMatrix = readRDS("functions/plotConfusionMatrix.rds")
  plotVip = readRDS("functions/plotVip.rds")
  
  
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
  
  #################################
  # Classification - parallelized
  #################################

  
  # choose number of cores - this code chooses half the cores available on your
  # machine
  numCores = floor(detectCores()/2)
  cluster = makeCluster(numCores)
  iterations = seq(1, numIterations)
  
  # this function will run PLSDA and save model objects into the specified directory
  # the console will log a list of NULL values, don't worry about this
  
  parLapply(cl = cluster, iterations, runPlsda, spectra = spectra,
            className = className, ncomp = ncomps, resampling = "down",
            include_age = includeAge, modelDirectory = '', saveModelObject = FALSE,
            cmDirectory = "", accuracyDirectory = accuracyDirectory)
  
  # get optimal number of components from downsampled model objects
  optimalNumComps = getNComps(accuracyDirectory, ncomps)
  
  # run plsda with upsampling with optimal number of components
  parLapply(cl = cluster, iterations, runPlsda, spectra = spectra,
            className = className, ncomp = optimalNumComps, resampling = "up",
            include_age = includeAge, modelDirectory = upSamplingDirectory, 
            cmDirectory = cmDirectory, saveModelObject = TRUE)
  
  ##############################################################################
  # Get data and plot confusion matrices and variable importance values
  ##############################################################################
  
  # get overall accuracy, mean confusion matrix, and standard deviation (2) 
  # confusion matrix
  
  data = getData(directory = cmDirectory, metricsDirectory = metricsDirectory,
                 className = className, includesAge = includeAge)
  
  matrixName = paste(paste(className, age, data$ncomps, 'comps', sep = '_'),
                     'jpeg', sep = '.')
  
  #plot confusion matrix as high resolution jpeg
  plotConfusionMatrix(data$cmMean, directory = cmFinalDirectory,
                      fileName = matrixName)
  
  #plot top 5 and bottom 5 variable importance values
  plotVip(modelDirectory = upSamplingDirectory, saveDirectory = vipDirectory,
          baseFileName = paste(className, 'vip', sep = '_'))
}

saveRDS(generatePLSData, 'functions/generatePLSData.rds')
