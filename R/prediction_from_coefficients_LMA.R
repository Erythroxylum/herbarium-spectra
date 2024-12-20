#'------------------------------------------------------------------------------
#' @title Predicting LMA using coefficients and validation datasets
#'------------------------------------------------------------------------------

#' @description A script to predict traits with validation data from coefficents
#' 
#' @return 1) *obs_pred.csv for plotting predicted and observed values
#' 2) *_performance.csv detailing statistics of model performance.

#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------

library(data.table)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("auxiliary/pred_coef_bands.R")
source("auxiliary/model_performance_average.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# SET the root_folder to read data and export results

# Antonio
root_path <- getwd()
results <- paste0(root_path, "/results")

#'------------------------------------------------------------------------------
#' @Read-Information_Reshape
#-------------------------------------------------------------------------------
# SET source
coef_source <- "Kothari" #Kothari or HUH
spectra_source <- "HUH" #Kothari or HUH
dataset <- "cwt" # ref, cwt, refnorm, cwtnorm, refnormcut, cwt

# SET bands of interest
bands <- seq(450, 2400, by = 5)
bands <- seq(1350, 2400, by = 5)
#bands <- seq(450, 1300, by = 5)

# remove sensor overlap region
#bands <- bands[!(bands >= 980 & bands <= 1000)]


# ------------------------------------------------------------------------------
###  Define out_path of results

out_path <- paste0(root_path, "/results/model_transfer/LMA_coef-", coef_source, "_spec-", spectra_source, "_", dataset, "_", min(bands), "-", max(bands))
dir.create(out_path, recursive = TRUE)

# ------------------------------------------------------------------------------
### Read coefficients

coef <- fread(paste0(root_path, "/results/", 
                     coef_source, "/", dataset, "/",coef_source, "_", min(bands), "-", max(bands), "/pls_", coef_source, "_coefficients.csv"))
colnames(coef) <- gsub("`", "", colnames(coef))

# ------------------------------------------------------------------------------
### Read spectra

if(spectra_source == "HUH") {
  
  # HUH
  if (dataset == "cwtnorm") {
    frame <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf561_cwt5nm_norm_450-2400.csv"))
  } else if (dataset == "cwt") {
    frame <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf561_cwt5nm_450-2400.csv"))
  } else if (dataset == "ref") {
    frame <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf561_ref5nm_450-2400.csv"))
  } else if (dataset == "refnorm") {
    frame <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf561_ref5nm_norm_450-2400.csv"))
  } else {
    stop("Invalid dataset specified for HUH.")
  }
  
  # drop rows with no trait data
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
  
} else if(spectra_source == "Kothari") {
  
  # Kothari
  if (dataset == "cwtnorm") {
    frame <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_cwt5nm_norm_450-2400.csv"))
  } else if (dataset == "refnorm") {
    frame <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_ref5nm_norm_450-2400.csv"))
  } else if (dataset == "cwt") {
    frame <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_cwt5nm_450-2400.csv"))
  } else if (dataset == "ref") {
    frame <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_ref5nm_450-2400.csv"))
  } else {
    stop("Invalid dataset specified for Kothari.")
  }
  
  # Drop rows with no trait data
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
#' @Predict-from-spectra
#-------------------------------------------------------------------------------

trait_predicted <- pred_coef(spectra, coef, bands)

#-------------------------------------------------------------------------------
#' @LMA_Performance
#-------------------------------------------------------------------------------

# Estimate performance
performance <- pbmclapply(X = 1:nrow(coef),
                          FUN = application_performance,
                          meta_split = meta,
                          trait_observed = frame$leafKg_m2,
                          predicted = trait_predicted,
                          mc.preschedule = TRUE, 
                          mc.set.seed = FALSE,
                          mc.cores = 2,
                          mc.cleanup = TRUE)

# Make it a frame
performance <- do.call(rbind, performance)
performance <- cbind(iteration = 1:nrow(performance), 
                     nsamp = nrow(frame),
                     performance)

# Export
fwrite(performance, paste0(out_path, "/coef-", coef_source, "_spec-", spectra_source, "_", dataset,
                           "_", min(bands), "-", max(bands), "_performance.csv"))

#-------------------------------------------------------------------------------
#' @LMA_Observed-predicted
#-------------------------------------------------------------------------------

# Get residuals
residuals <- frame$leafKg_m2 - trait_predicted

# Residuals statistics
residuals_mean <- rowMeans(residuals)
residuals_sd <- apply(residuals, 1, sd)
residuals_se <- residuals_sd / sqrt(nrow(coef))

# Predictions statistics
predicted_mean <- rowMeans(trait_predicted)
predicted_sd <- apply(trait_predicted, 1, sd)
predicted_se <- predicted_sd / sqrt(nrow(coef))

# Combine with meta information and summarize by `accession`
predicted_summary <- meta
predicted_summary$observed <- frame$leafKg_m2
predicted_summary$predicted_mean <- predicted_mean
predicted_summary$predicted_sd <- predicted_sd
predicted_summary$predicted_se <- predicted_se
predicted_summary$residuals_mean <- residuals_mean
predicted_summary$residuals_sd <- residuals_sd
predicted_summary$residuals_se <- residuals_se

# Export
fwrite(predicted_summary, paste0(out_path, "/coef-", coef_source, "_spec-", spectra_source, 
                                 "_", dataset, "_", min(bands), "-", max(bands), "_obs-pred.csv"))

#-------------------------------------------------------------------------------
#' @LMA_Visualization_basic
#-------------------------------------------------------------------------------

plot(predicted_summary$predicted_mean, predicted_summary$observed, main=paste0("coeff:",coef_source, ", spectra:", spectra_source, ", dataset:", dataset, "\nrange:", min(bands),"-",max(bands)))
abline(a = 0, b = 1, lty = "dotted", col = "grey")
abline(a = mean(performance$intercept), b = mean(performance$slope), col = "red")

pdf(paste0(out_path, "/From_", coef_source, "_to_", spectra_source, "_", dataset, "_", min(bands), "-", max(bands), "_obs-pred.pdf"), width = 8, height = 8)
plot(predicted_summary$predicted_mean, predicted_summary$observed, main=paste0("coeff:",coef_source, ", spectra:", spectra_source, ", dataset:", dataset, "\nrange:", min(bands),"-",max(bands)))
abline(a = 0, b = 1, lty = "dotted", col = "grey")
abline(a = mean(performance$intercept), b = mean(performance$slope), col = "red")
dev.off()

#-------------------------------------------------------------------------------
#' @Clear-environment
#-------------------------------------------------------------------------------

# remove non-function objects from envt
rm(list = ls()[!sapply(ls(), function(x) is.function(get(x)))])

#-------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------





