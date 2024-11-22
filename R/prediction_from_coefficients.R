#'------------------------------------------------------------------------------
#' @title Predicting leaf traits using coefficients
#'------------------------------------------------------------------------------

#' @description A script to predict traits from coefficents
#' 
#' @return several .csv files with according to the leaf trait of interest

#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------

library(data.table)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("auxiliary/pred_coef.R")
source("auxiliary/model_performance_average.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results

# Antonio
root_path <- getwd()
results <- paste0(root_path, "/results/ref")

#'------------------------------------------------------------------------------
#' @Read-Information_Reshape
#-------------------------------------------------------------------------------

# Define bands of interest
# bands <- seq(450, 2400, by = 5)
# bands <- seq(1350, 2400, by = 5)
# bands <- seq(450, 1300, by = 5)
bands <- seq(680, 900, by = 5)

# Define source
coef_source <- "HUH" #Kothari or HUH
spectra_source <- "Kothari" #Kothari or HUH

# Define out_path of results
out_path <- paste0(root_path, "/results/ref/model_transfer/From_", coef_source, "_to_", spectra_source, 
                   "_", min(bands), "-", max(bands))
dir.create(out_path, recursive = TRUE)

# ------------------------------------------------------------------------------
### Read coefficients

coef <- fread(paste0(root_path, "/results/ref/", 
                     coef_source, "/", coef_source, "_", min(bands), "-", max(bands),
                     "/pls_", coef_source, "_coefficients.csv"))
colnames(coef) <- gsub("`", "", colnames(coef))

# ------------------------------------------------------------------------------
### Read spectra

if(spectra_source == "HUH") {
  
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
  
} else if(spectra_source == "Kothari") {
  
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
#' @Predict-from-spectra
#-------------------------------------------------------------------------------

trait_predicted <- pred_coef(spectra, coef)

#-------------------------------------------------------------------------------
#' @Performance
#-------------------------------------------------------------------------------

# Estimate performance
performance <- pbmclapply(X = 1:nrow(coef),
                          FUN = application_performance,
                          meta_split = meta,
                          trait_observed = frame$leafKg_m2,
                          predicted = trait_predicted,
                          mc.preschedule = TRUE, 
                          mc.set.seed = FALSE,
                          mc.cores = 25,
                          mc.cleanup = TRUE)

# Make it a frame
performance <- do.call(rbind, performance)
performance <- cbind(iteration = 1:nrow(performance), 
                     nsamp = nrow(frame),
                     performance)

# Export
fwrite(performance, paste0(out_path, "/From_", coef_source, "_to_", spectra_source, 
                           "_", min(bands), "-", max(bands), "_performance.csv"))

#-------------------------------------------------------------------------------
#' @Observed-predicted
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
fwrite(predicted_summary, paste0(out_path, "/From_", coef_source, "_to_", spectra_source, 
                                 "_", min(bands), "-", max(bands), "_obs-pred.csv"))

#-------------------------------------------------------------------------------
#' @Visualization (basic)
#-------------------------------------------------------------------------------

plot(predicted_summary$predicted_mean, predicted_summary$observed)
abline(a = 0, b = 1, lty = "dotted", col = "grey")
abline(a = mean(performance$intercept), b = mean(performance$slope), col = "red")

pdf(paste0(out_path, "/obs-pred_transfer.pdf"), width = 8, height = 8)
plot(predicted_summary$predicted_mean, predicted_summary$observed)
abline(a = 0, b = 1, lty = "dotted", col = "grey")
abline(a = mean(performance$intercept), b = mean(performance$slope), col = "red")
dev.off()

#-------------------------------------------------------------------------------
#' @END
#-------------------------------------------------------------------------------