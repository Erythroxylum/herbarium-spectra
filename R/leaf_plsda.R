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

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------
# for DW MBP, setwd()
source("auxiliary/data_split.R")
source("auxiliary/data_segments.R")
source("auxiliary/folds.R")
source("auxiliary/model_tune_plsda.R")
source("auxiliary/model_build_classification.R")
source("auxiliary/model_performance.R")
source("auxiliary/pls_coefficients.R")
source("auxiliary/pls_vip.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results

root_path <- "C:/Users/jog4076/Downloads"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

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

species <- frame$species
species <- sub(" ", "_", species)

spectra <- frame[, .SD, .SDcols = 22:ncol(frame)]

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Get rows for training
split <- data_split(meta = meta)

# Export for record
saveRDS(split, paste0(root_path, "/classification_split.rds"))

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
                       mc.cores = 1) # If windows = 1

# Export for record
saveRDS(segments, paste0(root_path, "/classification_segments.rds"))

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
                               threads = 1) # If windows = 1, mac 2 same as 6

# Select the optimal based on the accuracy
plot(x = 1:ncomp_max,
     y = colMeans(opt_models[metric == "Accuracy", 4:ncol(opt_models)]),
     xlab = "Number of components",
     ylab = "PRESS")

# Manual selection
ncomp <- 31
abline(v = ncomp, col = "red")

# Export for record and figures
fwrite(opt_models, paste0(root_path, "/opt_comp_models.csv"))

#-------------------------------------------------------------------------------
#' @Final_model for prediction
#-------------------------------------------------------------------------------

models <- model_build_plsda(meta = meta,
                            split = split,
                            segments = segments,
                            species = species,
                            spectra = spectra,
                            ncomp = ncomp,
                            threads = 1) # If windows = 1


#-------------------------------------------------------------------------------
#' @Final_model for prediction
#-------------------------------------------------------------------------------

models_plsda <- model_build_plsda(meta = meta,
                            split = split,
                            segments = segments,
                            species = species,
                            spectra = spectra,
                            ncomp = ncomp,
                            threads = 1) # If windows = 1

models_lda <- model_build_lda(meta = meta,
                              split = split,
                              segments = segments,
                              species = species,
                              spectra = spectra,
                              threads = 1) # If windows = 1


#-------------------------------------------------------------------------------
#' @Performance-training
#-------------------------------------------------------------------------------

# This return the stats of the model performance and the predicted probabilities

performance_plsda_training <- model_performance_plsda(meta_split = meta[split,],
                                                species_split = species[split], 
                                                spectra_split = spectra[split, ],
                                                models = models_plsda,
                                                ncomp = ncomp,
                                                threads = 1)

performance_plsda_testing <- model_performance_plsda(meta_split = meta[!split,],
                                               species_split = species[!split], 
                                               spectra_split = spectra[!split, ],
                                               models = models_plsda,
                                               ncomp = ncomp,
                                               threads = 1)

performance_lda_training <- model_performance_lda(meta_split = meta[split,],
                                                  species_split = species[split], 
                                                  spectra_split = spectra[split, ],
                                                  models = models_lda,
                                                  threads = 1)

performance_lda_testing <- model_performance_lda(meta_split = meta[!split,],
                                                 species_split = species[!split], 
                                                 spectra_split = spectra[!split, ],
                                                 models = models_lda,
                                                 threads = 1)





setwd("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/herbarium_spectra/kothari_model/")

# Load necessary libraries
library(spectrolab)
library(caret)
library(dplyr)
library(ggplot2)

##### read data, convert to spectra object, save RDS
spec_df <- read.csv("../../fullDataHUH2024_sp25leaf636_noResample_400-2400.csv", header = T, check.names = F)

herb_spec_all <- as_spectra(spec_df, name_idx = 1, meta_idxs = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,17,18))

#saveRDS(herb_spec_all, file = "spectra/sp10HUH2024_sp25leaf636_noResample_400-2300.rds")

# optionally, read in spectra
#herb_spec_all = readRDS("spectra/specHUH2024_1nm_sp25leaf639-18oct2024-sp10.rds")

#####################################################
## Shorten species names to G. species.
meta(herb_spec_all)$SpShort <- factor(
  unlist(lapply(strsplit(as.character(meta(herb_spec_all)$species), split=" "),
                function(x) paste(substring(x[[1]], 1, 1), x[[2]], sep = ". ")))
)

## shorten to Gen. species
meta(herb_spec_all)$SpShort <- factor(
  unlist(lapply(strsplit(as.character(meta(herb_spec_all)$species), split = " "),
                function(x) paste(substring(x[[1]], 1, 3), x[[2]], sep = ". ")))
)

##################################################
## filter species, or not

# To run all species:
sp <- names(table(meta(herb_spec_all)$SpShort))

# Or select only 10 'common' species in Kothari dataset
#sp <- c("A. flexuosa", "A. rubrum", "A. saccharinum", "A. saccharum", "B. papyrifera", "B. populifolia", "F. grandifolia", "P. grandidentata", "P. tremuloides", "Q. rubra")
sp <- c("Ago. flexuosa", "Ace. rubrum", "Ace. saccharinum", "Ace. saccharum", "Bet. papyrifera", "Bet. populifolia", "Fag. grandifolia", "Pop. grandidentata", "Pop. tremuloides", "Que. rubra")

# Filter the herb_spec_all object to include only the spectra for the species listed in sp.
herb_spec <- herb_spec_all[which(meta(herb_spec_all)$SpShort %in% sp)]

# Remove unused levels from the SpShort factor in the herb_spec object.
meta(herb_spec)$SpShort <- droplevels(meta(herb_spec)$SpShort)

# Assign the filtered SpShort abbreviations from the metadata of herb_spec to a new variable herb_spec_class.
herb_spec_class <- meta(herb_spec)$SpShort

# Create training and testing datasets (60% for training)
train_idx <- createDataPartition(
  y = meta(herb_spec)$SpShort,
  p = .6,
  list = FALSE
)

herb_spec_train <- as.matrix(herb_spec[train_idx])
herb_spec_test <- as.matrix(herb_spec[-train_idx])
herb_spec_class_train <- herb_spec_class[train_idx]
herb_spec_class_test <- herb_spec_class[-train_idx]

# Set up cross-validation controls with down-sampling
ctrl1 <- trainControl(method = "repeatedcv",
                      repeats = 10,
                      number = 10,
                      summaryFunction = multiClassSummary,
                      sampling = "down")

# Set up cross-validation controls with up-sampling
ctrl2 <- trainControl(method = "repeatedcv",
                      repeats = 10,
                      number = 10,
                      summaryFunction = multiClassSummary,
                      sampling = "up")

# Train a PLS-DA model using down-sampling
plsFit_herb1 <- train(
  x = herb_spec_train,
  y = herb_spec_class_train,
  method = "pls",
  tuneLength = 50,   # Tune for up to 50 components
  trControl = ctrl1,
  probMethod = "softmax"
)

# Plot accuracy vs. number of components
# Extract the optimal number of components and accuracy
optimal_ncomp <- plsFit_herb1$bestTune$ncomp
optimal_accuracy <- plsFit_herb1$results$Accuracy[plsFit_herb1$results$ncomp == optimal_ncomp]

# Create the plot
acc_plot <- ggplot(plsFit_herb1$results, aes(x = ncomp, y = Accuracy)) +
  geom_line() +
  geom_point() +
  # Highlight the optimal point
  geom_point(data = plsFit_herb1$results[plsFit_herb1$results$ncomp == optimal_ncomp,],
             aes(x = ncomp, y = Accuracy),
             color = "red", size = 3) +
  # Add the label inside the plot
  annotate("text", x = optimal_ncomp, y = optimal_accuracy,
           label = paste(optimal_ncomp, "components,\naccuracy =", round(optimal_accuracy, 3)),
           hjust = .8, vjust = 1.5, size = 5, color = "red") +
  labs(title = "Cross-validation Accuracy for Different Components",
       x = "Number of Components",
       y = "Accuracy") +
  theme_minimal() +
  # Increase the size of axis labels
  theme(axis.title = element_text(size = 15),    # Axis title size
        axis.text = element_text(size = 15))     # Axis text size
# Save the plot
plsda_dn_acc <- "s10comp50/PLSDA_acc_dn_s10comp50.pdf"
ggsave(plsda_dn_acc, plot = acc_plot, width = 6, height = 6)

# Train another PLS-DA model using up-sampling and optimal ncomp from previous model
plsFit_herb2 <- train(
  x = herb_spec_train,
  y = herb_spec_class_train,
  method = "pls",
  tuneGrid = expand.grid(ncomp = plsFit_herb1$bestTune$ncomp),
  trControl = ctrl2,
  probMethod = "softmax"
  )

# Confusion matrix for the test set
herb_conmat <- confusionMatrix(
  predict(plsFit_herb2, herb_spec_test),
  herb_spec_class_test
  )

# Create and save the confusion matrix
saveRDS(herb_conmat, file = "confusion_matrix_sp10comp50.rds")

# Data manipulation for plotting the confusion matrix
pcm_d <- as.data.frame(herb_conmat$table)
pcm_d$Prediction <- unlist(lapply(as.character(pcm_d$Prediction), function(x){
  pred_split <- strsplit(x, split = " ")
  pred_paste <- paste(pred_split[[1]][1:2], collapse = " ")
  return(pred_paste)
}))

pcm_d$Reference <- unlist(lapply(as.character(pcm_d$Reference), function(x){
  pred_split <- strsplit(x, split = " ")
  pred_paste <- paste(pred_split[[1]][1:2], collapse = " ")
  return(pred_paste)
}))

pcm_d$RefSum <- unlist(lapply(pcm_d$Reference, function(x) sum(pcm_d$Freq[pcm_d$Reference == x])))
pcm_d$RefPer <- round(pcm_d$Freq / pcm_d$RefSum * 100, digits = 0)



###############################################################
### Plot

# Define the color palette from white to orange
cols <- colorRampPalette(c('white', '#fe9929'))

# Plot using ggplot, similar to corrplot style, with light grid and no legend
pcm_d_p <- ggplot(data = pcm_d[pcm_d$RefPer > 0,],
                  aes(x = Prediction , y = Reference, fill = RefPer)) +
  geom_tile(color = "white") +  # Add white border to tiles
  geom_text(aes(label = RefPer), color = '#542788', size = 6) +  # Purple for text
  theme_minimal() +  # Simplified theme
  scale_fill_gradientn(colours = cols(10), limits = c(0, 100)) +  # Gradient from white to orange
  guides(fill = "none") +  # Remove the legend
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic", size = 17),
    axis.text.y = element_text(face = "italic", size = 17),
    panel.grid.major = element_line(color = "gray90"),  # Light major grid lines
    panel.grid.minor = element_line(color = "gray95"),  # Light minor grid lines
    panel.grid.major.x = element_line(color = "gray90"),  # Add vertical major grid lines
    panel.grid.minor.x = element_line(color = "gray95"),  # Add vertical minor grid lines
    panel.border = element_blank()  # Remove border
  ) +
  ggtitle("PLS-DA of herbarium leaf spectra") +
  labs(y = "True identity") # add tag="A" for fig

# Save the plot
plsda_cm <- "s10comp50/PLSDA_cm_s10comp50.pdf"
ggsave(plsda_cm, plot = pcm_d_p, width = 7.5, height = 7.5)

##################################

# Save the final PLS-DA model object
saveRDS(plsFit_herb2, file = "s10comp50/plsda_upmodel_final_s10comp50.rds")

