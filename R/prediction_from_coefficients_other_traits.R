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
library(pbmcapply)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("auxiliary/pred_coef_bands.R")
source("auxiliary/model_performance_average.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results

root_path <- getwd()
results <- paste0(root_path, "/results")

#'------------------------------------------------------------------------------
#' @Read-Information_Reshape
#-------------------------------------------------------------------------------

# Define source
coef_source <- "Kothari" # Coefficients of trait model
spectra_source <- "HUH" # Spectra for new predictions
dataset <- "ref" # cwt or ref
trait_name <- "" # "" for lma, other C, Ca, car, cel, chlA, N, sol

# Define bands of interest
bands <- seq(450, 2400, by = 5) # full-range
#bands <- seq(1350, 2400, by = 5) # SWIRs
#bands <- seq(450, 1300, by = 5) # VNIR

# remove sensor overlap region
#bands <- bands[!(bands >= 980 & bands <= 1000)]

# Define out_path of results
out_path <- paste0(root_path, "/results/model_transfer/", "Coef-", coef_source, "_Spec-", spectra_source, "_", dataset, trait_name, "_", min(bands), "-", max(bands))
dir.create(out_path, recursive = TRUE)

# ------------------------------------------------------------------------------
### Read coefficients

coef <- fread(paste0(root_path, "/results/", 
                     coef_source, "/", dataset, "/",coef_source, "_", min(bands), "-", max(bands), "/pls_", coef_source, "_",dataset,"_",min(bands), "-", max(bands), "_coefficients.csv"))
colnames(coef) <- gsub("`", "", colnames(coef))

# ------------------------------------------------------------------------------
### Read spectra

if(spectra_source == "HUH") {
  
  # HUH
  if (dataset == "cwtnorm") {
    frame <- fread(paste0(root_path, "/data/DMWhiteHUHspec1_sp25leaf560_cwt5nm_norm_450-2400.csv"))
  } else if (dataset == "cwt") {
    frame <- fread(paste0(root_path, "/data/DMWhiteHUHspec1_sp25leaf560_cwt5nm_450-2400.csv"))
  } else if (dataset == "ref") {
    frame <- fread(paste0(root_path, "/data/DMWhiteHUHspec1_sp25leaf560_ref5nm_450-2400.csv"))
  } else if (dataset == "refnorm") {
    frame <- fread(paste0(root_path, "/data/DMWhiteHUHspec1_sp25leaf560_ref5nm_norm_450-2400.csv"))
  } else {
    stop("Invalid dataset specified for HUH.")
  }
  
  # drop rows with no trait data
  #frame <- frame[!is.na(leafKg_m2),]
  
  # HUH meta
  meta <- frame[, c("collector", "specimenIdentifier", "targetTissueClass", "targetTissueNumber", "measurementIndex", "scientificName",
                    "Genus", "Family", "Class", "Order",
                    "eventDate", "Age", "measurementFlags",
                    "tissueNotes", "hasGlue", "tissueDevelopmentalStage", "greenIndex", "growthForm", "doyOfCollection")]

  meta$sample <- 1:nrow(meta)
  
  # HUH traits
  traits <- frame[, c("leafKg_m2")]
  
  # HUH bands
  cbands <- as.character(bands)
  spectra <- frame[, ..cbands]
  
} else if(spectra_source == "Kothari") {
  
  # Kothari
  if (dataset == "cwtnorm") {
    frame <- fread(paste0(root_path, "/data/KothariPressedSpec_pressed_unavg_cwt5nm_norm_450-2400.csv"))
  } else if (dataset == "refnorm") {
    frame <- fread(paste0(root_path, "/data/KothariPressedSpec_pressed_unavg_ref5nm_norm_450-2400.csv"))
  } else if (dataset == "cwt") {
    frame <- fread(paste0(root_path, "/data/KothariPressedSpec_pressed_unavg_cwt5nm_450-2400.csv"))
  } else if (dataset == "ref") {
    frame <- fread(paste0(root_path, "/data/KothariPressedSpec_pressed_unavg_ref5nm_450-2400.csv"))
  } else {
    stop("Invalid dataset specified for Kothari.")
  }
  
  # Drop rows with no trait data
  #frame <- frame[!is.na(leafKg_m2),]
  
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
#' @Output_predictions
#-------------------------------------------------------------------------------

# Compute mean and standard deviation across all iterations for each sample
trait_stats <- trait_predicted[, .(
  predicted_trait_mean = rowMeans(.SD, na.rm = TRUE),
  predicted_trait_sd = apply(.SD, 1, sd, na.rm = TRUE)
)]

# Combine the computed statistics with metadata columns (1:22)
results_trait <- cbind(frame[, 1:22], trait_stats) #huh
#results_trait <- cbind(frame[, 1:59], trait_stats) #kot

print(results_trait)

# Optionally write to a CSV file
write.csv(results_trait, paste0(out_path, "/Coef-", coef_source, "_Spec-", spectra_source, "_", dataset, trait_name, "_", min(bands), "-", max(bands), "_predicted_trait.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
#' @Clear-environment
#-------------------------------------------------------------------------------

# remove non-function objects from envt
rm(list = ls()[!sapply(ls(), function(x) is.function(get(x)))])

#-------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------


