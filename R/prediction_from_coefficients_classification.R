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
library(parallel)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------
# for DW MBP, setwd()
source("auxiliary/confusion_matrices.R")
source("auxiliary/pred_coef.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results

#DW
root_path <- getwd()

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# SET source
coef_source <- "HUH" #Kot or HUH
spectra_source <- "Kot" #Kot or HUH
dataset <- "ref" # ref, cwt, refnorm, cwtnorm, refnormcut, cwt

# SET bands of interest
bands <- seq(450, 2400, by = 5)

# ------------------------------------------------------------------------------
###  Define out_path of results
out_path <- paste0(root_path, "/out_classify/model_transfer/LMA_coef-", coef_source, "_spec-", spectra_source, "_", dataset, "_", min(bands), "-", max(bands))

# Ensure the folder exists
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}


# ------------------------------------------------------------------------------
### Read spectra

if(spectra_source == "HUH") {
  
  # HUH
  if (dataset == "cwt") {
    frame <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf561_cwt5nm_450-2400.csv"))
  } else if (dataset == "ref") {
    frame <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf561_ref5nm_450-2400.csv"))
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
  if (dataset == "cwt") {
    frame <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_cwt5nm_450-2400.csv"))
  } else if (dataset == "ref") {
    frame <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_ref5nm_450-2400.csv"))
  } else {
    stop("Invalid dataset specified for Kothari.")
  }
  
  # Kothari meta
  meta <- frame[, c("name", "accession", "Species", "LatinGenus", "LatinSpecies", 
                    "Project", "Discoloration", "Stage", "growthForm", "Latitude", "Longitude",
                    "species", "sp10", "greenIndex")]
  
  #meta <- meta[sp10 == TRUE,]
  meta$sample <- 1:nrow(meta)
  
  # Kothari bands
  cbands <- as.character(bands)
  spectra <- frame[, ..cbands]
  
}

species_vector <- c("Quercus rubra", "Populus tremuloides", "Populus grandidentata", "Fagus grandifolia", "Betula populifolia", "Betula papyrifera", "Agonis flexuosa", "Acer saccharum", "Acer saccharinum", "Acer rubrum")

# Subset the data.table based on the species vector
frame <- frame[species %in% species_vector]

#-------------------------------------------------------------------------------
#' @Data_reshape
#-------------------------------------------------------------------------------

# Get files from meta data, traits, and spectra.

meta <- frame[, c("collector", "accession", "accession_leaf", "leaf", "scan","species",
                  "genus","family","class","order",
                  "ddmmyyScanned", "absoluteAge", "herbQuality",
                  "damage", "glue", "leafStage", "greenIndex")]

meta$sample <- 1:nrow(meta)

species <- frame$species

species <- sub(" ", "_", species)

spectra <- frame[, .SD, .SDcols = 23:ncol(frame)]


#-------------------------------------------------------------------------------
#' @Predict-from-spectra
#-------------------------------------------------------------------------------

## LDA
# Define classes from training data
classes <- c(
  'Acer_rubrum', 'Acer_saccharinum', 'Acer_saccharum', 
  'Agonis_flexuosa', 'Betula_papyrifera', 'Betula_populifolia', 
  'Fagus_grandifolia', 'Populus_grandidentata', 'Populus_tremuloides', 
  'Quercus_rubra'
)

# Read coefficients
coef <- fread(paste0(root_path, "/out_classify_MBP/coefficients_lda_sp10.csv"))

# Predict classes for new spectral data
predicted_classes <- predict_lda(spectra, coefficients, num_lds = 9, classes = classes)

# View predicted classes
print(predicted_classes)

# View predictions
print(predictions$probabilities)  # Probabilities for each class
print(predictions$classes)        # Final predicted classes

## PLSDA

classes <- c(
  'Acer_rubrum', 'Acer_saccharinum', 'Acer_saccharum', 
  'Agonis_flexuosa', 'Betula_papyrifera', 'Betula_populifolia', 
  'Fagus_grandifolia', 'Populus_grandidentata', 'Populus_tremuloides', 
  'Quercus_rubra'
)

# New spectral dataset
new_spectra <- fread("data/dataKothari_pressed_unavg_ref5nm_450-2400.csv")

# Predict classes for the new spectral dataset
predictions <- predict_plsda(new_spectra, coefficients, classes)

# View predictions
print(predictions$probabilities)  # Probabilities for each class
print(predictions$classes)        # Final predicted classes



#-------------------------------------------------------------------------------
#' @Performance
#-------------------------------------------------------------------------------

model_performance_transfer <- function(meta_split, species_split, spectra_split, models, ncomp, threads = 1) {
  # Predictions for all iterations
  predicted_results <- pred_coef(spectra_split, models, ncomp)
  
  # Initialize performance metrics
  accuracy_list <- list()
  precision_list <- list()
  balanced_accuracy_list <- list()
  
  # Iterate over predictions
  for (i in seq_along(predicted_results)) {
    predicted_classes <- predicted_results[[i]]$classes
    
    # Calculate confusion matrix
    cm <- table(Predicted = predicted_classes, Actual = species_split)
    
    # Calculate metrics
    accuracy <- sum(diag(cm)) / sum(cm)
    precision <- mean(diag(cm) / rowSums(cm))  # Per-class precision
    balanced_accuracy <- mean(diag(cm) / colSums(cm))  # Balanced accuracy
    
    # Append metrics
    accuracy_list[[i]] <- accuracy
    precision_list[[i]] <- precision
    balanced_accuracy_list[[i]] <- balanced_accuracy
  }
  
  # Aggregate metrics across iterations
  performance_summary <- list(
    mean_accuracy = mean(unlist(accuracy_list)),
    sd_accuracy = sd(unlist(accuracy_list)),
    mean_precision = mean(unlist(precision_list)),
    sd_precision = sd(unlist(precision_list)),
    mean_balanced_accuracy = mean(unlist(balanced_accuracy_list)),
    sd_balanced_accuracy = sd(unlist(balanced_accuracy_list))
  )
  
  return(performance_summary)
}

# This return the stats of the model performance and the predicted probabilities

performance_plsda <- model_performance_transfer(meta_split = meta,
                                                      species_split = species, 
                                                      spectra_split = spectra,
                                                      models = models_plsda,
                                                      ncomp = ncomp,
                                                      threads = 1)
#-------------------------------------------------------------------------------
#' @Confusion-Matrices
#-------------------------------------------------------------------------------

# Confusion Matrix Function
confusion_matrices_transfer <- function(meta_split, species_split, spectra_split, models, ncomp, threads = 1) {
  # Predictions for all iterations
  predicted_results <- pred_coef(spectra_split, models, ncomp)
  
  # Aggregate confusion matrices
  confusion_matrices <- lapply(predicted_results, function(prediction) {
    table(Predicted = prediction$classes, Actual = species_split)
  })
  
  # Combine confusion matrices into one
  combined_cm <- Reduce("+", confusion_matrices)
  
  return(combined_cm)
}

# Save Confusion Matrix to CSV
save_confusion_matrix <- function(confusion_matrix, output_path) {
  cm_df <- as.data.table(as.table(confusion_matrix))  # Convert to data.table for saving
  fwrite(cm_df, file = output_path)
}

CM_plsda <- confusion_matrices_transfer(meta_split = meta,
                                                 species_split = species, 
                                                 spectra_split = spectra,
                                                 models = models_plsda,
                                                 ncomp = ncomp,
                                                 threads = 2)

# Save Confusion Matrix to CSV
save_confusion_matrix(confusion_matrix, paste0(root_path, "/confusion_matrix_plsda.csv"))