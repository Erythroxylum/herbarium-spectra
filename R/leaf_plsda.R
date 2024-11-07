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
source("auxiliary/confusion_matrices_dw.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results

#JAG
root_path <- "C:/Users/jog4076/Downloads"
#DW
root_path <- getwd()

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

# Select the number of specimens per species to include in training
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
                               threads = 2) # If windows = 1, mac 2 same as 6

# Plot Accuracy
pdf("../herbarium_spectra_results/ncomp_accuracy.pdf", width = 5, height = 5) 
# Select the optimal based on the accuracy
plot(x = 1:ncomp_max,
     y = colMeans(opt_models[metric == "Accuracy", 4:ncol(opt_models)]),
     xlab = "Number of components",
     ylab = "Accuracy") # was PRESS
# Calculate mean accuracy for each component
mean_accuracy <- colMeans(opt_models[metric == "Accuracy", 4:ncol(opt_models)])
# Identify the component with the highest mean accuracy
ncomp <- which.max(mean_accuracy)
max_accuracy <- mean_accuracy[ncomp]
# Add a point and annotate the highest accuracy point
points(ncomp, max_accuracy, col = "red", pch = 19, cex = 1.2)
text(ncomp, max_accuracy, labels = paste("Max Acc.:", round(max_accuracy, 2)),
     pos = 1, offset = 2, col = "red")
text(ncomp, max_accuracy, labels = paste("N comp.:", ncomp),
     pos = 1, offset = 3, col = "red")
# Close the PDF device to save the file
dev.off()

# Manual selection
#ncomp <- 35
#abline(v = ncomp, col = "red")

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
                                  threads = 2) # If windows = 1

# Based on LDA
models_lda <- model_build_lda(meta = meta,
                              split = split,
                              segments = segments,
                              species = species,
                              spectra = spectra,
                              threads = 2) # If windows = 1


#-------------------------------------------------------------------------------
#' @Performance-training
#-------------------------------------------------------------------------------

# This return the stats of the model performance and the predicted probabilities

performance_plsda_training <- model_performance_plsda(meta_split = meta[split, ],
                                                      species_split = species[split], 
                                                      spectra_split = spectra[split, ],
                                                      models = models_plsda,
                                                      ncomp = ncomp,
                                                      threads = 2)

#generate inverse of numeric vector for species split
inverse_split <- setdiff(1:length(species), split)

performance_plsda_testing <- model_performance_plsda(meta_split = meta[!split, ],
                                                     species_split = species[inverse_split], 
                                                     spectra_split = spectra[!split, ],
                                                     models = models_plsda,
                                                     ncomp = ncomp,
                                                     threads = 1)

performance_lda_training <- model_performance_lda(meta_split = meta[split, ],
                                                  species_split = species[split], 
                                                  spectra_split = spectra[split, ],
                                                  models = models_lda,
                                                  threads = 2)

performance_lda_testing <- model_performance_lda(meta_split = meta[!split, ],
                                                  species_split = species[inverse_split], 
                                                  spectra_split = spectra[!split, ],
                                                  models = models_lda,
                                                  threads = 2)

# Export for record
saveRDS(performance_plsda_training, paste0(root_path, "/performance_plsda_training.rds"))
saveRDS(performance_plsda_testing, paste0(root_path, "/performance_plsda_testing.rds"))
saveRDS(performance_lda_training, paste0(root_path, "/performance_lda_training.rds"))
saveRDS(performance_lda_testing, paste0(root_path, "/performance_lda_testing.rds"))

#-------------------------------------------------------------------------------
#' @Confusion-Matrices
#-------------------------------------------------------------------------------

CM_plsda_training <- confusion_matrices_plsda_dw(meta_split = meta[split,],
                                             species_split = species[split], 
                                             spectra_split = spectra[split, ],
                                             models = models_plsda,
                                             ncomp = ncomp,
                                             threads = 1)

CM_plsda_testing <- confusion_matrices_plsda_dw(meta_split = meta[!split,],
                                            species_split = species[inverse_split], 
                                            spectra_split = spectra[!split, ],
                                            models = models_plsda,
                                            ncomp = ncomp,
                                            threads = 1)

CM_lda_training <- confusion_matrices_lda_dw(meta_split = meta[split,],
                                         species_split = species[split], 
                                         spectra_split = spectra[split, ],
                                         models = models_lda,
                                         threads = 1)

CM_lda_testing <- confusion_matrices_lda_dw(meta_split = meta[!split,],
                                        species_split = species[inverse_split], 
                                        spectra_split = spectra[!split, ],
                                        models = models_lda,
                                        threads = 1)

# Export
saveRDS(CM_plsda_training, paste0(root_path, "/CM_plsda_training.rds"))
saveRDS(CM_plsda_testing, paste0(root_path, "/CM_plsda_testing.rds"))
saveRDS(CM_lda_training, paste0(root_path, "/CM_lda_training.rds"))
saveRDS(CM_lda_testing, paste0(root_path, "/CM_lda_testing.rds"))

