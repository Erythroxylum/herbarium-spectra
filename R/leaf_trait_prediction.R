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
library(caret)
library(plsVarSel)
library(parallel)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("auxiliary/data_split_growthForm.R")
source("auxiliary/data_segments.R")
source("auxiliary/model_tune_sigma.R")
source("auxiliary/model_build_pls.R")
source("auxiliary/model_performance_average.R")
source("auxiliary/pls_coefficients.R")
source("auxiliary/pls_vip.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results
# Antonio
#root_path <- "C:/Users/jog4076/Downloads"
# DW
root_path <- "/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/herbarium-spectra"


#'------------------------------------------------------------------------------
#' @Read-Information_Reshape
#-------------------------------------------------------------------------------

# Select the file of interest

# HUH
frame <- fread("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/dataHUH2024_sp25leaf563_norm_5nm_400-2400.csv")
frame <- frame[!is.na(leafKg_m2),]
# HUH data
meta <- frame[, c("collector", "accession", "accession_leaf", "leaf", "scan", 
                  "species","genus","family","order","class",
                  "growthForm", "ddmmyyScanned","doy","absoluteAge", 
                  "herbQuality","damage", "glue", "leafStage", "greenIndex")]
meta$sample <- 1:nrow(meta)
traits <- frame[, c("leafKg_m2", "leafThickness")]
spectra <- frame[, .SD, .SDcols = 23:ncol(frame)]


# Kothari
frame <- fread("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/dataKothari_pressed_unavg_norm_5nm_400-2400.csv")
meta <- frame[, c("accession", "speciesAuthor", "species", "genus", "LatinEpithet", 
                  "sp10", "Project", "Discoloration", "Stage", "growthForm", 
                  "Latitude", "Longitude", "normalization_magnitude", "greenIndex")]
meta$sample <- 1:nrow(meta)
traits <- frame[, c("leafKg_m2", "LDMC", "EWT", "EWT_rehydrated", "N", "C", "NDF", 
                    "ADF", "ADL", "solubles", "hemicellulose", "cellulose", "lignin", 
                    "chlA", "chlB", "car", "Al", "Ca", "Cu", "Fe", "K", "Mg", "Mn", 
                    "Na", "P", "Zn", "N_area", "C_area", "solubles_area", 
                    "hemicellulose_area", "cellulose_area", "lignin_area", 
                    "chlA_area", "chlB_area", "car_area", "Al_area", "Ca_area", 
                    "Cu_area", "Fe_area", "K_area", "Mg_area", "Mn_area", "Na_area", 
                    "P_area", "Zn_area")]
spectra <- frame[, .SD, .SDcols = 61:ncol(frame)]


#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Get rows for training 
split <- data_split(meta = meta)

# Export for record
saveRDS(split, paste0(root_path, "/pls_huh5nm_split.rds"))
saveRDS(split, paste0(root_path, "/pls_kothari5nm_split.rds"))

#-------------------------------------------------------------------------------
#' @Segments
#-------------------------------------------------------------------------------

# Select a spectral measurement per specimen
iterations <- 10
segments <- pbmclapply(X = 1:iterations,
                       FUN = data_segments,
                       meta = meta,
                       split = split,
                       mc.set.seed = TRUE,
                       mc.cores = 1) # If windows = 1

# Export for record
saveRDS(segments, paste0(root_path, "/pls_huh5nm_segments.rds"))
saveRDS(segments, paste0(root_path, "/pls_kothari5nm_segments.rds"))

#-------------------------------------------------------------------------------
#' @Model_tune
#-------------------------------------------------------------------------------

# Select the optimal number of components for all the iterations

# Max number of comp to run, can't run over 24.
ncomp_max <- 20

# Models to evaluate the optimal
opt_models <- model_tune(meta = meta,
                         split = split,
                         segments = segments,
                         traits = traits[,1],
                         spectra = spectra,
                         ncomp_max = ncomp_max,
                         threads = 1) # If windows = 1

#select mode of Ncomp across models (from one-sigma method)
ncomp <- as.numeric(names(sort(table(opt_models$ncomp), decreasing = TRUE)[1]))

write(ncomp, "ncomp_onesigma-herbarium.txt")

# Calculate optimal number of components
n_components_PRESS <- which.min(colMeans(opt_models[metric == "PRESS" | estimate == "PRESS", 6:ncol(opt_models)]))
n_components_RMSEP <- which.min(colMeans(opt_models[metric == "RMSEP" | estimate == "RMSEP", 6:ncol(opt_models)]))
n_components_R2 <- which.max(colMeans(opt_models[metric == "R2" | estimate == "R2", 6:ncol(opt_models)]));

# Print all values in a single line
cat(sprintf("Optimal number of components based on max or min values - PRESS: %d, RMSEP: %d, R2: %d\n", 
            n_components_PRESS, n_components_RMSEP, n_components_R2), 
    file = "ncomp_otherMetrics-herbarium.txt")

# Plot optima
#PRESS
pdf("Ncomp_PRESS_plot.pdf", width = 5, height = 4)
plot(x = 1:(ncol(opt_models)-5),
     y = colMeans(opt_models[metric == "PRESS" | estimate == "PRESS", 6:ncol(opt_models)]),
     main = "Optimal number of components",
     xlab = "Number of components",
     ylab = "PRESS")
dev.off()

#R2
pdf("Ncomp_R2_plot.pdf", width = 5, height = 4)
plot(x = 1:(ncol(opt_models)-5),
     y = colMeans(opt_models[metric == "R2" | estimate == "R2", 6:ncol(opt_models)]),
     main = "Optimal number of components",
     xlab = "Number of components",
     ylab = "R2")
dev.off()

#RMSEP
pdf("Ncomp_RMSEP_plot.pdf", width = 5, height = 4)
plot(x = 1:(ncol(opt_models)-5),
     y = colMeans(opt_models[metric == "RMSEP" | estimate == "RMSEP", 6:ncol(opt_models)]),
     main = "Optimal number of components",
     xlab = "Number of components",
     ylab = "RMSEP")
dev.off()

# Export csv of statistics for record and figures
fwrite(opt_models, paste0(root_path, "/pls_huh5nm_opt_comp_models.csv"))
#fwrite(opt_models, paste0(root_path, "/pls_kothari5nm_opt_comp_models.csv"))

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

# Export models
saveRDS(models, paste0(root_path, "/pls_huh5nm_final_models.rds"))
#saveRDS(models, paste0(root_path, "/pls_kothari5nm_final_models.rds"))

#-------------------------------------------------------------------------------
#' @Performance-training
#-------------------------------------------------------------------------------

# Herbarium
training_performance <- model_performance(meta_split = meta[split,],
                                          traits_split = traits[split, 1], 
                                          spectra_split = spectra[split,],
                                          models = models,
                                          ncomp = ncomp,
                                          threads = 1) # If windows = 1

fwrite(training_performance$performance, paste0(root_path, "/pls_huh5nm_training_performance.csv"))
fwrite(training_performance$predicted, paste0(root_path, "/pls_huh5nm_training_obs-pred.csv"))
fwrite(training_performance$predicted_by_accession, paste0(root_path, "/pls_huh5nm_training_obs-pred-avg.csv"))

# Pressed
training_performance_kothari <- model_performance_kothari(meta_split = meta[split,],
                                                          traits_split = traits[split, 1], 
                                                          spectra_split = spectra[split,],
                                                          models = models,
                                                          ncomp = ncomp,
                                                          threads = 1) # If windows = 1

# Export for record and figures
fwrite(training_performance_kothari$performance, paste0(root_path, "/pls_kothari5nm_training_performance.csv"))
fwrite(training_performance_kothari$predicted, paste0(root_path, "/pls_kothari5nm_training_obs-pred.csv"))
fwrite(training_performance_kothari$predicted_by_accession, paste0(root_path, "/pls_kothari5nm_training_obs-pred-avg.csv"))

#-------------------------------------------------------------------------------
#' @Performance-testing
#-------------------------------------------------------------------------------

# Herbarium
testing_performance <- model_performance(meta_split = meta[!split,],
                                         traits_split = traits[!split, 1], 
                                         spectra_split = spectra[!split,],
                                         models = models,
                                         ncomp = ncomp,
                                         threads = 1) # If windows = 1

# write
fwrite(testing_performance$performance, paste0(root_path, "/pls_huh5nm_testing_performance.csv"))
fwrite(testing_performance$predicted, paste0(root_path, "/pls_huh5nm_testing_obs-pred.csv"))
fwrite(testing_performance$predicted_by_accession, paste0(root_path, "/pls_huh5nm_testing_obs-pred-avg.csv"))

#Pressed
testing_performance_kothari <- model_performance_kothari(meta_split = meta[!split,],
                                         traits_split = traits[!split, 1], 
                                         spectra_split = spectra[!split,],
                                         models = models,
                                         ncomp = ncomp,
                                         threads = 1) # If windows = 1

fwrite(testing_performance_kothari$performance, paste0(root_path, "/pls_kothari5nm_testing_performance.csv"))
fwrite(testing_performance_kothari$predicted, paste0(root_path, "/pls_kothari5nm_testing_obs-pred.csv"))
fwrite(testing_performance_kothari$predicted_by_accession, paste0(root_path, "/pls_kothari5nm_testing_obs-pred-avg.csv"))

#-------------------------------------------------------------------------------
#' @Model_coefficients
#-------------------------------------------------------------------------------

coefficients <- pls_coefficients(models = models,
                                 ncomp = ncomp)

fwrite(coefficients, paste0(root_path, "/pls_huh5nm_coefficients.csv"))
#fwrite(coefficients, paste0(root_path, "/pls_kothari5nm_coefficients.csv"))

#-------------------------------------------------------------------------------
#' @Model_VIP
#-------------------------------------------------------------------------------

vip <- pls_vip(models = models,
               ncomp = ncomp)

fwrite(vip, paste0(root_path, "/pls_huh5nm_vip.csv"))
#fwrite(vip, paste0(root_path, "/pls_kothari5nm_vip.csv"))

#-------------------------------------------------------------------------------
#' @Pressed-vs-herbarium-performance
#-------------------------------------------------------------------------------

# pressMod_herbSpec: pressed models on herbarium spectra

# press models
pres_models <- readRDS("pls_kothari5nm_final_models.rds")

# load herbarium spectra and metadata 
frame <- fread("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/dataHUH2024_sp25leaf563_norm_5nm_400-2400.csv")
frame <- frame[!is.na(leafKg_m2),]
# HUH data
meta <- frame[, c("collector", "accession", "accession_leaf", "leaf", "scan", 
                  "species","genus","family","order","class",
                  "growthForm", "ddmmyyScanned","doy","absoluteAge", 
                  "herbQuality","damage", "glue", "leafStage", "greenIndex")]
meta$sample <- 1:nrow(meta)
traits <- frame[, c("leafKg_m2", "leafThickness")]
spectra <- frame[, .SD, .SDcols = 23:ncol(frame)]

# Run predictions, need ncomp from above
presMod_herbSpec <- model_performance(meta_split = meta,
                                         traits_split = traits[,1], 
                                         spectra_split = spectra,
                                         models = pres_models,
                                         ncomp = 9,
                                         threads = 1) # If windows = 1

fwrite(presMod_herbSpec$performance, paste0(root_path, "/pls_presMod_herbSpec_performance.csv"))
fwrite(presMod_herbSpec$predicted, paste0(root_path, "/pls_presMod_herbSpec_obs-pred.csv"))
fwrite(presMod_herbSpec$predicted_by_accession, paste0(root_path, "/pls_presMod_herbSpec_obs-pred-avg.csv"))


#-------------------------------------------------------------------------------
# herbMod_presSpec: herbarium models on pressed spectra

# herb models
herb_models <- readRDS("pls_huh5nm_final_models.rds")

# press spec
frame <- fread("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/dataKothari_pressed_unavg_norm_5nm_400-2400.csv")
meta <- frame[, c("accession", "speciesAuthor", "species", "genus", "LatinEpithet", 
                  "sp10", "Project", "Discoloration", "Stage", "growthForm", 
                  "Latitude", "Longitude", "normalization_magnitude", "greenIndex")]
meta$sample <- 1:nrow(meta)
traits <- frame[, c("leafKg_m2", "LDMC", "EWT", "EWT_rehydrated", "N", "C", "NDF", 
                    "ADF", "ADL", "solubles", "hemicellulose", "cellulose", "lignin", 
                    "chlA", "chlB", "car", "Al", "Ca", "Cu", "Fe", "K", "Mg", "Mn", 
                    "Na", "P", "Zn", "N_area", "C_area", "solubles_area", 
                    "hemicellulose_area", "cellulose_area", "lignin_area", 
                    "chlA_area", "chlB_area", "car_area", "Al_area", "Ca_area", 
                    "Cu_area", "Fe_area", "K_area", "Mg_area", "Mn_area", "Na_area", 
                    "P_area", "Zn_area")]
spectra <- frame[, .SD, .SDcols = 61:ncol(frame)]

## Run predictions, need ncomp from above
herbMod_presSpec <- model_performance_kothari(meta_split = meta,
                                                traits_split = traits[,1], 
                                                spectra_split = spectra,
                                                models = herb_models,
                                                ncomp = 10,
                                                threads = 1) # If windows = 1

# Write output
fwrite(herbMod_presSpec$performance, paste0(root_path, "/pls_herbMod_presSpec_performance.csv"))
fwrite(herbMod_presSpec$predicted, paste0(root_path, "/pls_herbMod_presSpec_obs-pred.csv"))
fwrite(herbMod_presSpec$predicted_by_accession, paste0(root_path, "/pls_herbMod_presSpec_obs-pred-avg.csv"))


#-------------------------------------------------------------------------------
#' @END
#-------------------------------------------------------------------------------
