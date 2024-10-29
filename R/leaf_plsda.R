#'------------------------------------------------------------------------------
#' @title Classification of herbarium specimens using leaf spectra
#'------------------------------------------------------------------------------

#' @description A script to train and evaluate PLS-DA classification of leaf
#' spectra to taxa
#'
#' @return several .csv files

#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------

library(data.table)
library(pls)
library(caret)
library(plsVarSel)
library(parallel)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------
# for DW MBP, setwd()
source("auxiliary/data_split.R")
source("auxiliary/data_segments.R")
source("auxiliary/folds.R")
source("auxiliary/model_tune_plsda.R")
source("auxiliary/model_build_classification.R")
source("auxiliary/model_performance_classification.R")
source("auxiliary/confusionMatrices.R")
#source("auxiliary/pls_coefficients.R")
#source("auxiliary/pls_vip.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results

root_path <- "C:/Users/jog4076/Downloads"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Select the file of interest
frame <- fread(paste0(root_path,
                      "/fullDataHUH2024_sp25leaf636_noResample_400-2300.csv")) 
frame <- frame[!is.na(leafKg_m2),]

#-------------------------------------------------------------------------------
#' @Data_reshape
#-------------------------------------------------------------------------------

# Get files from meta data, traits, and spectra.

meta <- frame[, c("collector", "accession", "accession_leaf", "leaf", "scan","species",
                  "ddmmyyScanned", "absoluteAge", "herbQuality",
                  "damage", "glue", "leafStage", "greenIndex")]

meta$sample <- 1:nrow(meta)

species <- frame$species
species <- sub(" ", "_", species)

spectra <- frame[, .SD, .SDcols = 22:ncol(frame)]

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Get rows for training
split <- data_split(meta = meta)

# Export for record
saveRDS(split, paste0(root_path, "/classification_split.rds"))

#-------------------------------------------------------------------------------
#' @Segments
#-------------------------------------------------------------------------------

# Select a spectral measurement per specimen
iterations <- 3 # 1000
segments <- pbmclapply(X = 1:iterations,
                       FUN = data_segments,
                       meta = meta,
                       split = split,
                       mc.set.seed = TRUE,
                       mc.cores = 1) # If windows = 1

# Export for record
saveRDS(segments, paste0(root_path, "/classification_segments.rds"))

#-------------------------------------------------------------------------------
#' @Model_tune
#-------------------------------------------------------------------------------

# Select the optimal number of components for all the iterations

# Max number of comp to run
ncomp_max <- 50

# Models to evaluate the optimal
opt_models <- model_tune_plsda(meta = meta,
                               split = split,
                               segments = segments,
                               species = species,
                               spectra = spectra,
                               ncomp_max = ncomp_max,
                               threads = 1) # If windows = 1, mac 2 same as 6

# Select the optimal based on the accuracy
plot(x = 1:ncomp_max,
     y = colMeans(opt_models[metric == "Accuracy", 4:ncol(opt_models)]),
     xlab = "Number of components",
     ylab = "PRESS")

# Manual selection
ncomp <- 39
abline(v = ncomp, col = "red")

# Export for record and figures
fwrite(opt_models, paste0(root_path, "/opt_comp_models.csv"))

#-------------------------------------------------------------------------------
#' @Final_model for prediction
#-------------------------------------------------------------------------------

# Based on PLSDA
models_plsda <- model_build_plsda(meta = meta,
                                  split = split,
                                  segments = segments,
                                  species = species,
                                  spectra = spectra,
                                  ncomp = ncomp,
                                  threads = 1) # If windows = 1

# Based on LDA
models_lda <- model_build_lda(meta = meta,
                              split = split,
                              segments = segments,
                              species = species,
                              spectra = spectra,
                              threads = 1) # If windows = 1


#-------------------------------------------------------------------------------
#' @Performance-training
#-------------------------------------------------------------------------------

# This return the stats of the model performance and the predicted probabilities

performance_plsda_training <- model_performance_plsda(meta_split = meta[split,],
                                                      species_split = species[split], 
                                                      spectra_split = spectra[split, ],
                                                      models = models_plsda,
                                                      ncomp = ncomp,
                                                      threads = 1)

performance_plsda_testing <- model_performance_plsda(meta_split = meta[!split,],
                                                     species_split = species[!split], 
                                                     spectra_split = spectra[!split, ],
                                                     models = models_plsda,
                                                     ncomp = ncomp,
                                                     threads = 1)

performance_lda_training <- model_performance_lda(meta_split = meta[split,],
                                                  species_split = species[split], 
                                                  spectra_split = spectra[split, ],
                                                  models = models_lda,
                                                  threads = 1)

performance_lda_testing <- model_performance_lda(meta_split = meta[!split,],
                                                 species_split = species[!split], 
                                                 spectra_split = spectra[!split, ],
                                                 models = models_lda,
                                                 threads = 1)

#-------------------------------------------------------------------------------
#' @Confusion-Matrices
#-------------------------------------------------------------------------------

CM_plsda_training <- confusionMatrices_plsda(meta_split = meta[split,],
                                             species_split = species[split], 
                                             spectra_split = spectra[split, ],
                                             models = models_plsda,
                                             ncomp = ncomp,
                                             threads = 1)

CM_plsda_testing <- confusionMatrices_plsda(meta_split = meta[!split,],
                                            species_split = species[!split], 
                                            spectra_split = spectra[!split, ],
                                            models = models_plsda,
                                            ncomp = ncomp,
                                            threads = 1)

CM_lda_training <- confusionMatrices_lda(meta_split = meta[split,],
                                         species_split = species[split], 
                                         spectra_split = spectra[split, ],
                                         models = models_lda,
                                         threads = 1)

CM_lda_testing <- confusionMatrices_lda(meta_split = meta[!split,],
                                        species_split = species[!split], 
                                        spectra_split = spectra[!split, ],
                                        models = models_lda,
                                        threads = 1)
