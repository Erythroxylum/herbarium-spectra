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
# for DW MBP, setwd()
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

## Set the path to export results, ensuring it is relative to root_path
output_path <- file.path(root_path, "out_classify")
# Ensure the output folder exists
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

## Set path for split and segments files
split_path <- file.path(root_path, "out_classify")
# Ensure the split folder exists
if (!dir.exists(split_path)) {
  dir.create(split_path, recursive = TRUE)
}

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Select the file of interest
frame <- fread("DMWhiteHUHspec1_sp25leaf560_ref5nm_450-2400.csv")

#-------------------------------------------------------------------------------
#' @Data_reshape
#-------------------------------------------------------------------------------

# Get files from meta data, traits, and spectra.

meta <- frame[, c("collector", "specimenIdentifier", "targetClass", "targetTissueNumber", "measurementIndex", "scientificName",
                  "Genus", "Family", "Class", "Order",
                  "eventDate", "Age", "measurementFlags",
                  "tissueNotes", "hasGlue", "tissueDevelopmentalStage", "greenIndex")]

# create sample index column
meta$sample <- 1:nrow(meta)

# set taxonomic level for classification: frame$scientificName or frame$genus
taxon <- frame$scientificName

# remove space
taxon <- sub(" ", "_", taxon)

# define spectra
spectra <- frame[, .SD, .SDcols = 23:ncol(frame)]

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Select the number of specimens per species to include in training
#split <- data_split(meta = meta)

# Export for record, reload for subsequent runs
#saveRDS(split, paste0(split_path, "/classification_split.rds"))

# Reload for comparability
split <- readRDS(paste0(split_path, "/classification_split.rds"))

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
                       mc.cores = 5) # If windows = 1

# Export for record
#saveRDS(segments, paste0(split_path, "/classification_segments.rds"))

# Reload for comparability
segments <- readRDS(paste0(split_path, "/classification_segments.rds"))

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
                               species = taxon,
                               spectra = spectra,
                               ncomp_max = ncomp_max,
                               threads = 20) # If windows = 1, mac 2 same as 6


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
ncomp_sd <- as.integer(which(mean_accuracy >= accuracy_threshold)[1])

#ncomp <- ncomp_sd
ncomp <- max_accuracy_component

# Print ncomp summary
cat(
  "The highest mean accuracy:", max_accuracy, "at component", max_accuracy_component, "\n",
  "Optimal component accuracy within one SD:", accuracy_threshold, "at component", ncomp_sd, "\n", "used:", ncomp, 
  file = paste0(output_path,"/ncomp_plsda.txt"), append = TRUE
)

# Export for record and figures
fwrite(opt_models, paste0(output_path, "/opt_comp_models_plsda.csv"))

#-------------------------------------------------------------------------------
#' @Final_model for prediction
#-------------------------------------------------------------------------------

# Based on PLSDA
models_plsda <- model_build_plsda(meta = meta,
                                  split = split,
                                  segments = segments,
                                  species = taxon,
                                  spectra = spectra,
                                  ncomp = ncomp,
                                  threads = 20) # If windows = 1

# save models
#saveRDS(models_plsda, output_path, "classification_plsda_models.rds")

#-------------------------------------------------------------------------------
#' @Performance-training
#-------------------------------------------------------------------------------

# This return the stats of the model performance and the predicted probabilities

performance_plsda_training <- model_performance_plsda(meta_split = meta[split, ],
                                                      species_split = taxon[split], 
                                                      spectra_split = spectra[split, ],
                                                      models = models_plsda,
                                                      ncomp = ncomp,
                                                      threads = 20)

#generate inverse of numeric vector for species split
inverse_split <- setdiff(1:length(taxon), split)

performance_plsda_testing <- model_performance_plsda(meta_split = meta[!split, ],
                                                     species_split = taxon[inverse_split], 
                                                     spectra_split = spectra[!split, ],
                                                     models = models_plsda,
                                                     ncomp = ncomp,
                                                     threads = 20)
# Export for record
saveRDS(performance_plsda_training, paste0(output_path, "/performance_plsda_training.rds"))
saveRDS(performance_plsda_testing, paste0(output_path, "/performance_plsda_testing.rds"))

#-------------------------------------------------------------------------------
#' @Model_coefficients
#-------------------------------------------------------------------------------

coefficients <- plsda_coefficients(models = models_plsda,
                                   ncomp = ncomp)

saveRDS(coefficients, paste0(output_path, "/coefficients_plsda.rds"))

#-------------------------------------------------------------------------------
#' @Model_VIP
#-------------------------------------------------------------------------------

vip_plsda <- plsda_vip(models = models_plsda)

fwrite(vip_plsda, paste0(output_path, "/varImp_plsda.csv"))

#-------------------------------------------------------------------------------
#' @Confusion-Matrices
#-------------------------------------------------------------------------------

CM_plsda_training <- confusion_matrices_plsda_dw(meta_split = meta[split,],
                                             species_split = taxon[split], 
                                             spectra_split = spectra[split, ],
                                             models = models_plsda,
                                             ncomp = ncomp,
                                             threads = 2)

CM_plsda_testing <- confusion_matrices_plsda_dw(meta_split = meta[!split,],
                                            species_split = taxon[inverse_split], 
                                            spectra_split = spectra[!split, ],
                                            models = models_plsda,
                                            ncomp = ncomp,
                                            threads = 2)
# Export
saveRDS(CM_plsda_training, paste0(output_path, "/CM_plsda_training.rds"))
saveRDS(CM_plsda_testing, paste0(output_path, "/CM_plsda_testing.rds"))

