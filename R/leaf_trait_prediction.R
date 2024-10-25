#'------------------------------------------------------------------------------
#' @title Predicting leaf traits from herbarium specimens
#'------------------------------------------------------------------------------

#' @description A script to train, evaluate, and predict leaf traits from 
#' herbarium specimens.
#' 
#' @return several .csv files with according to the leaf trait of interest

#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------

library(data.table)
library(pls)
library(plsVarSel)
library(parallel)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("auxiliary/data_split.R")
source("auxiliary/data_segments.R")
source("auxiliary/model_tune.R")
source("auxiliary/model_build.R")
source("auxiliary/model_performance.R")
source("auxiliary/pls_coefficients.R")
source("auxiliary/pls_vip.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results
# Antonio
root_path <- "C:/Users/jog4076/Downloads"
# DW
root_path <- "/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/herbarium-spectra"


#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# download from Google Drive
library(googledrive)
# drive_auth() authenticate for first use

# Replace with your actual file ID (https://drive.google.com/file/d/1RciaJCLduoDEO8bQWcGtZWaVxD9HA4V8/view?usp=drive_link)
file_id <- "1RciaJCLduoDEO8bQWcGtZWaVxD9HA4V8"  # Example file ID

# Specify the full path to save file
path <- "downloaded_file.csv"  # Adjust for your OS

# Download the file
drive_download(as_id(file_id), path = path, overwrite = TRUE)

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

traits <- frame[, c("leafKg_m2", "leafThickness")]

rank <- frame[, c("class", "order", "family", "genus", "spShort")]

spectra <- frame[, .SD, .SDcols = 22:ncol(frame)]


#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Get rows for training 
split <- data_split(meta = meta,
                    p = 0.6)

# Export for record
saveRDS(split, paste0(root_path, "/split.rds"))

#-------------------------------------------------------------------------------
#' @Segments
#-------------------------------------------------------------------------------

# Select a spectral measurement per specimen
iterations <- 1000
segments <- pbmclapply(X = 1:iterations,
                       FUN = data_segments,
                       meta = meta,
                       split = split,
                       mc.set.seed = TRUE,
                       mc.cores = 1) # If windows = 1

# Export for record
saveRDS(segments, paste0(root_path, "/segments.rds"))

#-------------------------------------------------------------------------------
#' @Model_tune
#-------------------------------------------------------------------------------

# Select the optimal number of components for all the iterations

# Max numeber of comp to run
ncomp_max <- 30

# Models to evaluate the optimal
opt_models <- model_tune(meta = meta,
                         split = split,
                         segments = segments,
                         traits = traits[,1],
                         spectra = spectra,
                         ncomp_max = ncomp_max,
                         threads = 1) # If windows = 1

# Select the optimal based on the CV behavior or PRESS values
#PRESS
plot(x = 1:ncomp_max,
     y = colMeans(opt_models[metric == "PRESS" | estimate == "PRESS", 6:ncol(opt_models)]),
     main = "Optimal number of components",
     xlab = "Number of components",
     ylab = "PRESS")

#R2
plot(x = 1:ncomp_max,
                    y = colMeans(opt_models[metric == "R2" | estimate == "R2", 6:ncol(opt_models)]),
                    main = "Optimal number of components",
                    xlab = "Number of components",
                    ylab = "R2")

#RMSEP
pdf("Ncomp_RMSEP_plot.pdf", width = 5, height = 4)
plot(x = 1:ncomp_max,
     y = colMeans(opt_models[metric == "RMSEP" | estimate == "RMSEP", 6:ncol(opt_models)]),
     main = "Optimal number of components",
     xlab = "Number of components",
     ylab = "RMSEP",
     abline(v = ncomp, col = "red"))
dev.off()

 # Calculate optimal number of components
n_components_PRESS <- which.min(colMeans(opt_models[metric == "PRESS" | estimate == "PRESS", 6:ncol(opt_models)]))
n_components_RMSEP <- which.min(colMeans(opt_models[metric == "RMSEP" | estimate == "RMSEP", 6:ncol(opt_models)]))
n_components_R2 <- which.max(colMeans(opt_models[metric == "R2" | estimate == "R2", 6:ncol(opt_models)]));

# Print all values in a single line
cat(sprintf("Optimal number of components - PRESS: %d, RMSEP: %d, R2: %d\n", n_components_PRESS, n_components_RMSEP, n_components_R2))

# Manual selection
ncomp <- 14
abline(v = ncomp, col = "red")

# Export for record and figures
fwrite(opt_models, paste0(root_path, "/opt_comp_models.csv"))

#-------------------------------------------------------------------------------
#' @Final_model for prediction
#-------------------------------------------------------------------------------

models <- model_build(meta = meta,
                      split = split,
                      segments = segments,
                      traits = traits[,1],
                      spectra = spectra,
                      ncomp = ncomp,
                      threads = 1) # If windows = 1

#-------------------------------------------------------------------------------
#' @Performance-training
#-------------------------------------------------------------------------------

training_performance <- model_performance(meta_split = meta[split,],
                                          traits_split = traits[split, 1], 
                                          spectra_split = spectra[split,],
                                          models = models,
                                          ncomp = ncomp,
                                          threads = 1) # If windows = 1

# Export for record and figures
fwrite(training_performance$performance, paste0(root_path, "/training_performance.csv"))
fwrite(training_performance$predicted, paste0(root_path, "/training_obs-pred.csv"))

#-------------------------------------------------------------------------------
#' @Performance-testing
#-------------------------------------------------------------------------------

testing_performance <- model_performance(meta_split = meta[!split,],
                                         traits_split = traits[!split, 1], 
                                         spectra_split = spectra[!split,],
                                         models = models,
                                         ncomp = ncomp,
                                         threads = 1) # If windows = 1

fwrite(testing_performance$performance, paste0(root_path, "/testing_performance.csv"))
fwrite(testing_performance$predicted, paste0(root_path, "/testing_obs-pred.csv"))

#-------------------------------------------------------------------------------
#' @Model_coefficients
#-------------------------------------------------------------------------------

coefficients <- pls_coefficients(models = models,
                                 ncomp = ncomp)

fwrite(coefficients, paste0(root_path, "/pls_coefficients.csv"))

#-------------------------------------------------------------------------------
#' @Model_VIP
#-------------------------------------------------------------------------------

vip <- pls_vip(models = models,
               ncomp = ncomp)

fwrite(vip, paste0(root_path, "/pls_vip.csv"))

#-------------------------------------------------------------------------------
#' @END
#-------------------------------------------------------------------------------
