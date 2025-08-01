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
library(plsVarSel)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("auxiliary/data_split_balanced.R")
source("auxiliary/data_segments.R")
source("auxiliary/folds.R")
source("auxiliary/model_tune_plsda.R")
source("auxiliary/model_build_classification.R")
source("auxiliary/model_performance_classification.R")
source("auxiliary/confusion_matrices.R")
source("auxiliary/pls_coefficients.R")
source("auxiliary/pls_vip.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_path to read data and export results
root_path <- getwd()

# Set the path to export results, ensuring it is relative to root_path
output_path <- file.path(root_path, "out_classify")
# Ensure the output folder exists
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Set path for split and segments files
split_path <- file.path(root_path, "out_classify")
# Ensure the split folder exists
if (!dir.exists(split_path)) {
  dir.create(split_path, recursive = TRUE)
}

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Select the file of interest
frame <- fread(paste0(root_path, "/data/DMWhiteHUHspec1_sp25leaf560_ref5nm_450-2400.csv"))

# remove space in scientificName
frame$scientificName <- sub(" ", "_", frame$scientificName)

#-------------------------------------------------------------------------------
#' @Data_reshape
#-------------------------------------------------------------------------------

# Get files from meta data, traits, and spectra.

meta <- frame[, c("collector", "specimenIdentifier", "targetTissueClass", "targetTissueNumber", "measurementIndex", "scientificName",
                  "Genus", "Family", "Class", "Order",
                  "eventDate", "Age", "measurementFlags",
                  "tissueNotes", "hasGlue", "tissueDevelopmentalStage", "greenIndex", "growthForm")]

# create sample index column or assign from idx_analysis
meta$sample <- 1:nrow(meta)

# Define bands of interest
bands <- seq(450, 2400, by = 5)
cbands <- as.character(bands)
# define spectra
spectra <- frame[, ..cbands]

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Select the number of specimens per taxon to include in training
#split <- data_split(meta = meta)

# reload from split created by plsda, for comparability
split <- readRDS(paste0(split_path, "/classification_split.rds"))

#-------------------------------------------------------------------------------
#' @Segments
#-------------------------------------------------------------------------------

# Select a spectral measurement per specimen
iterations <- 3 # 1000

segments <- pbmclapply(X = 1:iterations,
                       FUN = data_segments_classification,
                       meta = meta,
                       split = split,
                       mc.set.seed = TRUE,
                       mc.cores = 1) # If windows = 1

# Export for record
#saveRDS(segments, paste0(split_path, "/classification_segments.rds"))

# Reload from PLSDA segments, for comparability
#segments <- readRDS(paste0(split_path, "/classification_segments.rds"))


#-------------------------------------------------------------------------------
#' @Final_model for prediction
#-------------------------------------------------------------------------------

# Based on LDA
models_lda <- model_build_lda(meta = meta,
                              split = split,
                              segments = segments,
                              rank = scientificName,
                              spectra = spectra,
                              threads = 2) # If windows = 1

# save models (heavy)
#saveRDS(models_lda, output_path, "classification_lda_models.rds")

#-------------------------------------------------------------------------------
#' @Performance-training
#-------------------------------------------------------------------------------

# This returns the stats of the model performance and the predicted probabilities
# Change meta$scientificName to other rank if needed.

#generate inverse of numeric vector for species split
inverse_split <- setdiff(1:length(meta$scientificName), split)

performance_lda_training <- model_performance_lda(meta_split = meta[split, ],
                                                  species_split = meta$scientificName[split],
                                                  spectra_split = spectra[split, ],
                                                  models = models_lda,
                                                  threads = 2)

performance_lda_testing <- model_performance_lda(meta_split = meta[!split, ],
                                                  species_split = meta$scientificName[inverse_split], 
                                                  spectra_split = spectra[!split, ],
                                                  models = models_lda,
                                                  threads = 10)

# Export for record
saveRDS(performance_lda_training, paste0(output_path, "/performance_lda_training.rds"))
saveRDS(performance_lda_testing, paste0(output_path, "/performance_lda_testing.rds"))

#-------------------------------------------------------------------------------
#' @Model_coefficients
#-------------------------------------------------------------------------------

coefficients <- lda_coefficients(models = models_lda)

fwrite(coefficients, paste0(output_path, "/coefficients_lda.csv"))

#-------------------------------------------------------------------------------
#' @Model_VIP
#-------------------------------------------------------------------------------

vip_plsda <- lda_vip(models = models_lda)

fwrite(vip_plsda, paste0(output_path, "/varImp_lda.csv"))

#-------------------------------------------------------------------------------
#' @Confusion-Matrices
#-------------------------------------------------------------------------------

CM_lda_training <- confusion_matrices_lda(meta_split = meta[split,],
                                         species_split = meta$scientificName[split], 
                                         spectra_split = spectra[split, ],
                                         models = models_lda,
                                         threads = 2)

CM_lda_testing <- confusion_matrices_lda(meta_split = meta[!split,],
                                        species_split = meta$scientificName[inverse_split], 
                                        spectra_split = spectra[!split, ],
                                        models = models_lda,
                                        threads = 2)

# Export
saveRDS(CM_lda_training, paste0(output_path, "/CM_lda_training.rds"))
saveRDS(CM_lda_testing, paste0(output_path, "/CM_lda_testing.rds"))


#-------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------
