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
library(googledrive)

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
root_path <- "/Users/dawsonwhite/Library/Mobile\ Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/herbarium-spectra"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Select the file of interest
frame <- fread("../fullDataHUH2024_sp25leaf636_noResample_400-2300.csv") 
frame <- frame[!is.na(leafKg_m2),]

#-------------------------------------------------------------------------------
#' @Data_reshape
#-------------------------------------------------------------------------------

# Get files from meta data, traits, and spectra.

meta <- frame[, c("collector", "accession", "accession_leaf", "leaf", "scan",
                  "species","genus","family","class","order","ddmmyyScanned", 
                  "absoluteAge", "herbQuality","damage", "glue", "leafStage", 
                  "greenIndex")]

meta$sample <- 1:nrow(meta)

species <- frame$species
species <- sub(" ", "_", species)
genus <- frame$genus

spectra <- frame[, .SD, .SDcols = 22:ncol(frame)]

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Select the number of specimens per species to include in training
split <- data_split(meta = meta)

# Export for record
saveRDS(split, "../herbarium_spectra_results/leaf_plsda_it20/classification_split.rds")

#-------------------------------------------------------------------------------
#' @Segments
#-------------------------------------------------------------------------------

# Select a spectral measurement per specimen
iterations <- 20 # 1000

segments <- pbmclapply(X = 1:iterations,
                       FUN = data_segments,
                       meta = meta,
                       split = split,
                       mc.set.seed = TRUE,
                       mc.cores = 1) # If windows = 1

# Export for record
saveRDS(segments, "../herbarium_spectra_results/leaf_plsda_it20/classification_segments.rds")

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
                               threads = 8) # If windows = 1, mac 2 same as 6

# Filter for accuracy and accuracy SD metrics
accuracy_data <- opt_models[metric == "Accuracy", 4:ncol(opt_models)]
accuracy_sd_data <- opt_models[metric == "AccuracySD", 4:ncol(opt_models)]

# Calculate the mean accuracy and mean accuracy SD across iterations
mean_accuracy <- colMeans(accuracy_data, na.rm = TRUE)
mean_accuracy_sd <- colMeans(accuracy_sd_data, na.rm = TRUE)

# Find the component with the highest mean accuracy
max_accuracy <- max(mean_accuracy, na.rm = TRUE)
max_accuracy_component <- which.max(mean_accuracy)

# Calculate the threshold: highest mean accuracy minus one standard deviation
accuracy_threshold <- max_accuracy - mean_accuracy_sd[max_accuracy_component]

# Identify the lowest component within one SD of the highest accuracy
ncomp <- as.integer(which(mean_accuracy >= accuracy_threshold)[1])

# or assign ncomp as the highest accuracy
ncomp <- as.integer(which.max(mean_accuracy))

# Print results
cat("The highest mean accuracy:", max_accuracy, "at component", max_accuracy_component, "\n")
cat("Optimal component accuracy within one SD:", accuracy_threshold, "at component", ncomp, "\n")

# Export for record and figures
fwrite(opt_models, "../herbarium_spectra_results/leaf_plsda_it20/opt_comp_models.csv")

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
                                  threads = 8) # If windows = 1

# Based on LDA
models_lda <- model_build_lda(meta = meta,
                              split = split,
                              segments = segments,
                              species = species,
                              spectra = spectra,
                              threads = 8) # If windows = 1


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

