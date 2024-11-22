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
root_path <- getwd()

#'------------------------------------------------------------------------------
#' @Read-Information_Reshape
#-------------------------------------------------------------------------------

# Define source
source <- "Kothari" #Kothari or HUH

# Define bands of interest
# bands <- seq(450, 2400, by = 5)
# bands <- seq(1350, 2400, by = 5)
# bands <- seq(450, 1300, by = 5)
bands <- seq(680, 900, by = 5)

# Define out_path of results
out_path <- paste0(root_path, "/results/", source, "/", source, "_", min(bands), "-", max(bands))
dir.create(out_path, recursive = TRUE)

if(source == "HUH") {
  
  # HUH
  # frame <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf563_cwt_450-2400.csv")) #CWT
  frame <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf563_ref_400-2400.csv")) #ref
  frame <- frame[!is.na(leafKg_m2),]
  
  # HUH meta
  meta <- frame[, c("collector", "accession", "accession_leaf", "leaf", "scan", 
                    "species","genus","family","order","class",
                    "growthForm", "ddmmyyScanned","doy","absoluteAge", 
                    "herbQuality","damage", "glue", "leafStage", "greenIndex")]
  meta$sample <- 1:nrow(meta)
  
  # HUH traits
  traits <- frame[, c("leafKg_m2", "leafThickness")]
  
  # HUH bands
  cbands <- as.character(bands)
  spectra <- frame[, ..cbands]
  
} else if(source == "Kothari") {
  
  # Kothari
  # frame <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_cwt_450-2400.csv")) #CWT
  frame <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_ref_400-2400.csv")) #ref
  frame <- frame[!is.na(leafKg_m2),]
  
  # Kothari meta
  meta <- frame[, c("name", "accession", "Species", "LatinGenus", "LatinSpecies", 
                    "Project", "Discoloration", "Stage", "growthForm", "Latitude", "Longitude",
                    "species", "sp10", "greenIndex")]
  #meta <- meta[sp10 == TRUE,]
  meta$sample <- 1:nrow(meta)
  
  # Kothari traits
  traits <- frame[, c("leafKg_m2", "LDMC", "EWT", "EWT_rehydrated", "N", "C", "NDF", 
                      "ADF", "ADL", "solubles", "hemicellulose", "cellulose", "lignin", 
                      "chlA", "chlB", "car", "Al", "Ca", "Cu", "Fe", "K", "Mg", "Mn", 
                      "Na", "P", "Zn", "N_area", "C_area", "solubles_area", 
                      "hemicellulose_area", "cellulose_area", "lignin_area", 
                      "chlA_area", "chlB_area", "car_area", "Al_area", "Ca_area", 
                      "Cu_area", "Fe_area", "K_area", "Mg_area", "Mn_area", "Na_area", 
                      "P_area", "Zn_area")]
  
  # Kothari bands
  cbands <- as.character(bands)
  spectra <- frame[, ..cbands]
  
}

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

#' # Get rows for training
#' split <- data_split(meta = meta, p = 0.7)
#' 
#' # Export for record
#' saveRDS(split, paste0(out_path, "/pls_", source, "_split.rds"))

split <- readRDS(paste0(root_path, "/results/", source, "/pls_", source, "_split.rds"))

#' #-------------------------------------------------------------------------------
#' #' @Segments
#' #-------------------------------------------------------------------------------
#' 
#' # Select a spectral measurement per specimen
#' iterations <- 100
#' segments <- pbmclapply(X = 1:iterations,
#'                        FUN = data_segments,
#'                        meta = meta,
#'                        split = split,
#'                        mc.set.seed = TRUE,
#'                        mc.cores = 25) # If windows = 1
#' 
#' # Export for record
#' saveRDS(segments, paste0(out_path, "/pls_", source, "_segments.rds"))

segments <- readRDS(paste0(root_path, "/results/", source, "/pls_", source, "_segments.rds"))

#-------------------------------------------------------------------------------
#' @Model_tune
#-------------------------------------------------------------------------------

# Select the optimal number of components for all the iterations
ncomp_max <- 30

# Models to evaluate the optimal
opt_models <- model_tune(meta = meta,
                         split = split,
                         segments = segments,
                         traits = traits[,1],
                         spectra = spectra,
                         ncomp_max = ncomp_max,
                         threads = 25) # If windows = 1

# Plot optima
#PRESS
pdf(paste0(out_path, "/Ncomp_PRESS_plot.pdf"), width = 5, height = 4)
plot(x = 1:30,
     y = colMeans(opt_models[metric == "PRESS" | estimate == "PRESS", 7:ncol(opt_models)]),
     main = "Optimal number of components",
     xlab = "Number of components",
     ylab = "PRESS")
dev.off()

#R2
pdf(paste0(out_path, "/Ncomp_R2_plot.pdf"), width = 5, height = 4)
plot(x = 1:30,
     y = colMeans(opt_models[metric == "R2" | estimate == "R2", 7:ncol(opt_models)]),
     main = "Optimal number of components",
     xlab = "Number of components",
     ylab = "R2")
dev.off()

#RMSEP
pdf(paste0(out_path,"/Ncomp_RMSEP_plot.pdf"), width = 5, height = 4)
plot(x = 1:30,
     y = colMeans(opt_models[metric == "RMSEP" | estimate == "RMSEP", 7:ncol(opt_models)]),
     main = "Optimal number of components",
     xlab = "Number of components",
     ylab = "RMSEP")
dev.off()

# Manually select ncomp
ncomp <- 14

# Export csv of statistics for record and figures
fwrite(opt_models, paste0(out_path, "/pls_", source, "_opt_comp_models.csv"))
fwrite(data.table(ncomp = ncomp), paste0(out_path, "/", source, "_opt_comp.csv"))

#-------------------------------------------------------------------------------
#' @Final_model for prediction
#-------------------------------------------------------------------------------

models <- model_build(meta = meta,
                      split = split,
                      segments = segments,
                      traits = traits[,1],
                      spectra = spectra,
                      ncomp = ncomp,
                      threads = 25) # If windows = 1

# Export models
#saveRDS(models, paste0(out_path, "/pls_", source, "_final_models.rds"))

#-------------------------------------------------------------------------------
#' @Performance-training
#-------------------------------------------------------------------------------

training_performance <- model_performance(meta_split = meta[split,],
                                          traits_split = traits[split, 1], 
                                          spectra_split = spectra[split,],
                                          models = models,
                                          ncomp = ncomp,
                                          threads = 25) # If windows = 1

# Export
fwrite(training_performance$performance, paste0(out_path, "/pls_", source, "_training_performance.csv"))
fwrite(training_performance$predicted, paste0(out_path, "/pls_", source, "_training_obs-pred.csv"))

#-------------------------------------------------------------------------------
#' @Performance-testing
#-------------------------------------------------------------------------------

testing_performance <- model_performance(meta_split = meta[!split,],
                                         traits_split = traits[!split, 1], 
                                         spectra_split = spectra[!split,],
                                         models = models,
                                         ncomp = ncomp,
                                         threads = 25) # If windows = 1

# Export
fwrite(testing_performance$performance, paste0(out_path, "/pls_", source, "_testing_performance.csv"))
fwrite(testing_performance$predicted, paste0(out_path, "/pls_", source, "_testing_obs-pred.csv"))

#-------------------------------------------------------------------------------
#' @Model_coefficients
#-------------------------------------------------------------------------------

coefficients <- pls_coefficients(models = models,
                                 ncomp = ncomp)

fwrite(coefficients, paste0(out_path, "/pls_", source, "_coefficients.csv"))

#-------------------------------------------------------------------------------
#' @Model_VIP
#-------------------------------------------------------------------------------

vip <- pls_vip(models = models,
               ncomp = ncomp)

fwrite(vip, paste0(root_path, "/pls_", source, "_vip.csv"))

#-------------------------------------------------------------------------------
#' @END
#-------------------------------------------------------------------------------
