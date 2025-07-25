#'------------------------------------------------------------------------------
#' @title Regressions of classification probabilities against herbarium quality 
#' predictors
#'------------------------------------------------------------------------------

#' @description Script to classify spectra from PLS-DA coefficients, and run 
#' analyses testing herbarium specimen metadata against classification results
#' 
 
#' @return Plots
#' 

#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------

library(spectrolab)
library(caret)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(patchwork)
library(ape)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("auxiliary/predict_from_coefficients_classification.R")
source("auxiliary/confusion_matrices.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results
root_path <- getwd()
output_path <- paste0(root_path, "/Figures_Tables/")
# Ensure the output folder exists
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Read data
frame <- fread(paste0(root_path, 
                      "/data/DMWhiteHUHspec1_sp25leaf560_ref5nm_450-2400.csv"))

# remove space in scientificName
frame$scientificName <- sub(" ", "_", frame$scientificName)

#-------------------------------------------------------------------------------
#' @Data_reshape  
#-------------------------------------------------------------------------------

# Get files from meta data, traits, and spectra.

meta <- frame[, c("collector", "specimenIdentifier", "targetTissueClass", "targetTissueNumber", "measurementIndex", "scientificName",
                  "Genus", "Family", "Class", "Order",
                  "eventDate", "Age", "measurementFlags",
                  "tissueNotes", "hasGlue", "tissueDevelopmentalStage", "greenIndex", "growthForm", "doyOfCollection")]

# Convert IHerbSpec metadata fields to NPH nomenclature for leafDamage and specimenQuality
meta[, leafDamage := fcase(
  tissueNotes == "No visible biotic or abiotic damage to any leaves on herbarium sheet.", "none",
  tissueNotes == "Physical biotic or abiotic damage visible on some leaves on the specimen but no damage on the measured leaf.", "minor",
  tissueNotes == "Damage visible on measured leaves, but no damage is present in the  measured target area.", "medium",
  tissueNotes == "Biotic or abiotic damage is visible in the measured target area.", "major",
  default = NA_character_
)]
meta[, specimenQuality := gsub("Preservation\\|?.*", "", measurementFlags)]

traits <- frame[, c("leafKg_m2")]

meta$sample <- 1:nrow(meta)

species <- frame$scientificName

# Define bands of interest
bands <- seq(450, 2400, by = 5)
cbands <- as.character(bands)
# define spectra
spectra <- frame[, ..cbands]


#'------------------------------------------------------------------------------
#' @Predict_from_coefficients
#-------------------------------------------------------------------------------

# Read PLSDA coefficients
coef <- readRDS(paste0(root_path, "/out_classify/coefficients_plsda.rds"))

# Predict probabilities
predictions_table <- predict_plsda(
  coefficients_list = coef, 
  spectra = spectra, 
  meta = meta,
  traits = traits
  )

# write (but will add NTD and PD below)
#fwrite(predictions, paste0(output_path, "predictions_table_plsda25spp.csv"))


#'------------------------------------------------------------------------------
#' @Add-Phylogenetic-Distance
#-------------------------------------------------------------------------------

# Load the phylogenetic tree
phylo <- read.tree(paste0(root_path, "/phylogram_pd_TimeTree5.tre"))

# Compute the pairwise matrix of phylogenetic distances
distance_matrix <- cophenetic.phylo(phylo)

# Convert the distance matrix to a data.table
distance_matrix_dt <- as.data.table(as.table(distance_matrix), keep.rownames = "species1")

# Rename columns for clarity
setnames(distance_matrix_dt, c("V1", "V2", "N"), c("species1", "species2", "PD_value"))

# Ensure species names match the predictions data
distance_matrix_dt[, species1 := gsub(" ", "_", species1)]
distance_matrix_dt[, species2 := gsub(" ", "_", species2)]

# Calculate NTD for each true_class
ntd <- distance_matrix_dt[species1 != species2, .(
  ntd = min(PD_value, na.rm = TRUE)
), by = species1]

# Rename `species1` to `true_class` for merging compatibility
setnames(ntd, "species1", "true_class")

# Remove any duplicate or hidden attributes on 'predicted_class'
predictions_table <- as.data.table(predictions_table)
predictions_table[, predicted_class := as.character(predicted_class)]
predictions_table <- predictions_table[, !duplicated(names(predictions_table)), with = FALSE]

# Merge NTD into predictions
predictions_table_ntd <- merge(predictions_table, ntd, by = "true_class", all.x = TRUE)

# Add the phylogenetic distance to the predicted_class
predictions_table_ntd_pd <- merge(
  predictions_table_ntd, 
  distance_matrix_dt,
  by.x = c("true_class", "predicted_class"),
  by.y = c("species1", "species2"),
  all.x = TRUE
)

# reset pred and true as factors
predictions_table_ntd_pd$true_class <- as.factor(predictions_table_ntd_pd$true_class)
predictions_table_ntd_pd$predicted_class <- as.factor(predictions_table_ntd_pd$predicted_class)

# Rename the distance column for clarity
setnames(predictions_table_ntd_pd, "PD_value", "predicted_class_PD")

# Write
fwrite(predictions_table_ntd_pd, paste0(output_path, "predictions_table_plsda25spp.csv"))


#'------------------------------------------------------------------------------
#' @Subset_DataFrame
#-------------------------------------------------------------------------------

# Generate data frame with correct classifications
correct_classifications <- predictions_table_ntd_pd %>%
  filter(as.character(true_class) == as.character(predicted_class)) %>%
  mutate(dataset = "Correct Predictions")

# Generate data frame with incorrect classifications
incorrect_classifications <- predictions_table_ntd_pd %>%
  filter(as.character(true_class) != as.character(predicted_class)) %>%
  mutate(dataset = "Incorrect Predictions")

# Add a label for the full dataset
classProbs <- predictions_table_ntd %>%
  mutate(dataset = "All")


#'------------------------------------------------------------------------------
#' @Plots-Prob_vs_predictors-Fig7
#-------------------------------------------------------------------------------

#Correct classification bar plots
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr) # for ggarrange function
library(dunn.test)
library(ggsignif)

# Pick datasets
cor_inc_data <- bind_rows(correct_classifications, incorrect_classifications)

# Set y-axis limits for consistency
ymin = 0.03
ymax = 0.14

p1 <- ggplot(cor_inc_data, aes(x = specimenQuality, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  labs(x = "Specimen Quality", y = "Predicted Probabilities") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF", "#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Plot for glue
p2 <- ggplot(cor_inc_data, aes(x = hasGlue, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  labs(x = "Glue", y = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF", "#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Plot for damage - convert from tissueNotes
cor_inc_data[, leafDamage := fcase(
  tissueNotes == "No visible biotic or abiotic damage to any leaves on herbarium sheet.", "none",
  tissueNotes == "Physical biotic or abiotic damage visible on some leaves on the specimen but no damage on the measured leaf.", "minor",
  tissueNotes == "Damage visible on measured leaves, but no damage is present in the  measured target area.", "medium",
  tissueNotes == "Biotic or abiotic damage is visible in the measured target area.", "major",
  default = NA_character_
)]

p3 <- ggplot(cor_inc_data, aes(x = leafDamage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  labs(x = "Damage Status", y = "Predicted Probabilities") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF", "#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Plot for leafStage
p4 <- ggplot(cor_inc_data, aes(x = tissueDevelopmentalStage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  labs(x = "Leaf Stage", y = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF", "#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Arrange all plots in a grid
plots4 <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# Save the combined plot to a single png
ggsave("Figures_Tables/Fig7-PLSDA-regressions_correct-incorrect.png", plot = plots4, width = 8.5, height = 7)


######
## significance testing

#####
# specimenQuality
# Filter data for Correct Predictions
correct_data <- cor_inc_data[dataset == "Correct Predictions"]

# Plot for Specimen Quality within Correct Predictions
p1c <- ggplot(correct_data, aes(x = specimenQuality, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("good", "medium"), c("good", "poor"), c("medium", "poor")),  # Specify comparisons
    map_signif_level = TRUE,
    test = "t.test"
  ) +
  labs(x = "Specimen Quality", y = "Predicted Probabilities", 
       title = "Specimen Quality: Correct Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Filter data for Incorrect Predictions
incorrect_data <- cor_inc_data[dataset == "Incorrect Predictions"]

# Plot for Specimen Quality within Incorrect Predictions
p1i <- ggplot(incorrect_data, aes(x = specimenQuality, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("good", "medium"), c("good", "poor"), c("medium", "poor")),  # Specify comparisons
    map_signif_level = TRUE,
    test = "t.test"
  ) +
  labs(x = "Specimen Quality", y = "Predicted Probabilities", 
       title = "Specimen Quality: Incorrect Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax))

######
# Glue

# Filter data for Correct Predictions
correct_glue <- cor_inc_data[dataset == "Correct Predictions"]

# Filter data for Incorrect Predictions
incorrect_glue <- cor_inc_data[dataset == "Incorrect Predictions"]

# Plot for Glue (Correct Predictions)
p2_correct <- ggplot(correct_glue, aes(x = hasGlue, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("TRUE", "FALSE")),  # Replace with actual levels of glue
    map_signif_level = TRUE,
    test = "t.test"
  ) +
  labs(x = "Glue", y = "Predicted Probabilities", title = "Glue: Correct Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, by = 0.02))

# Plot for Glue (Incorrect Predictions)
p2_incorrect <- ggplot(incorrect_glue, aes(x = hasGlue, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("TRUE", "FALSE")),  # Replace with actual levels of glue
    map_signif_level = TRUE,
    test = "t.test"
  ) +
  labs(x = "Glue", y = "Predicted Probabilities", title = "Glue: Incorrect Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, by = 0.02))

#####
# leafDamage

# Filter data for Correct Predictions
correct_damage <- cor_inc_data[dataset == "Correct Predictions"]

# Filter data for Incorrect Predictions
incorrect_damage <- cor_inc_data[dataset == "Incorrect Predictions"]

# Plot for leafDamage (Correct Predictions)
p3_correct <- ggplot(correct_damage, aes(x = leafDamage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("major", "medium"), c("major","minor"),c("major","none"),c("medium","minor"),c("medium","none"),c("minor", "none")),  # Replace with actual levels of leafDamage
    map_signif_level = TRUE,
    test = "t.test"
  ) +
  labs(x = "Damage Status", y = "Predicted Probabilities", title = "Damage: Correct Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, by = 0.02))

p3_correct <- ggplot(correct_damage, aes(x = leafDamage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("major", "medium"), c("major", "minor"), c("major", "none"),
                       c("medium", "minor"), c("medium", "none"), c("minor", "none")),  # Levels of leafDamage
    map_signif_level = TRUE,
    test = "t.test",
    y_position = c(0.122, 0.126, 0.132, 0.112, 0.136, 0.122)  # Adjust y-positions to prevent overlap
  ) +
  labs(x = "Damage Status", y = "Predicted Probabilities", title = "Damage: Correct Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, by = 0.02))


# Plot for Damage (Incorrect Predictions)
p3_incorrect <- ggplot(incorrect_damage, aes(x = leafDamage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("major", "medium"), c("major", "minor"), c("major", "none"),
                       c("medium", "minor"), c("medium", "none"), c("minor", "none")),  # Levels of damage
    map_signif_level = TRUE,
    test = "t.test",
    y_position = c(0.122, 0.126, 0.132, 0.112, 0.136, 0.122)  # Adjust y-positions to prevent overlap
  ) +
  labs(x = "Damage Status", y = "Predicted Probabilities", title = "Damage: Incorrect Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, by = 0.02))


#####
# Leaf Stage

# Filter data for Correct Predictions
correct_leaf_stage <- cor_inc_data[dataset == "Correct Predictions"]

# Filter data for Incorrect Predictions
incorrect_leaf_stage <- cor_inc_data[dataset == "Incorrect Predictions"]

# Plot for Leaf Stage (Correct Predictions)
p4_correct <- ggplot(correct_leaf_stage, aes(x = leafStage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("mature", "young")),  # Replace with actual levels of leafStage
    map_signif_level = TRUE,
    test = "t.test"
  ) +
  labs(x = "Leaf Stage", y = "Predicted Probabilities", title = "Leaf Stage: Correct Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, by = 0.02))

# Plot for Leaf Stage (Incorrect Predictions)
p4_incorrect <- ggplot(incorrect_leaf_stage, aes(x = leafStage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("mature", "young")),  # Replace with actual levels of leafStage
    map_signif_level = TRUE,
    test = "t.test"
  ) +
  labs(x = "Leaf Stage", y = "Predicted Probabilities", title = "Leaf Stage: Incorrect Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, by = 0.02))


###### 1. Test for specimenQuality
kruskal_test_specimenQuality <- kruskal.test(prob_predicted ~ specimenQuality, data = classProbs)
print(kruskal_test_specimenQuality)

# Pairwise comparisons using Dunn's test
dunn_specimenQuality <- dunn.test(classProbs$prob_predicted, classProbs$specimenQuality, method = "bh")
dunn_specimenQuality
# [1] "good - medium" "good - poor"   "medium - poor"
# p 4.535799e-01 4.527629e-06 4.197685e-06

######  2. Test for glue
kruskal_test_glue <- kruskal.test(prob_predicted ~ hasGlue, data = correct_classifications)
print(kruskal_test_glue)
# Kruskal-Wallis chi-squared = 3.7555, df = 1, p-value = 0.05263

#######  3. Test for damage
kruskal_test_damage <- kruskal.test(prob_predicted ~ leafDamage, data = correct_classifications)
print(kruskal_test_damage)
#Kruskal-Wallis chi-squared = 15.65, df = 3, p-value = 0.001338

# Pairwise comparisons
dunn_damage <- dunn.test(correct_classifications$prob_predicted, correct_classifications$leafDamage, method = "bh")
print(dunn_damage$P.adjusted)
#$P.adjusted
#[1] 0.1505601144 0.2545721306 0.1375234946 0.1282708638 0.2344388191 0.0005490144
#$comparisons
#[1] "major - medium" "major - minor"  "medium - minor" "major - none"   "medium - none" 
#[6] "minor - none"  

#######  4. Test for leafStage
kruskal_test_leafStage <- kruskal.test(prob_predicted ~ tissueDevelopmentalStage, data = correct_classifications)
print(kruskal_test_leafStage)
# Kruskal-Wallis chi-squared = 0.0793, df = 1, p-value = 0.78

dev.off()



#'------------------------------------------------------------------------------
#' @Linear_Regressions_vs_numeric_predictors-Fig8
#-------------------------------------------------------------------------------

# Combine datasets for plotting
combined_data <- bind_rows(classProbs, correct_classifications, incorrect_classifications)

# Ensure dataset values align with color keys
colors <- c("Correct Predictions" = "#21908CFF", "Incorrect Predictions" = "#FDE725FF", "All" = "#440154FF")
regression_colors <- c("Correct Predictions" = "#116059FF", "Incorrect Predictions" = "#D4E735FF", "All" = "#440154FF")

# Plot Absolute Age vs Classification Probability
agelm <- ggplot(combined_data, aes(x = Age, y = prob_predicted, color = dataset)) +
  geom_point(size=0.9) +
  geom_smooth(method = "lm", se = FALSE, aes(color = dataset), linetype = "solid") +
  scale_color_manual(values = colors) +
  labs(title = NULL,
       x = "Age (years)",
       y = "Classification Probability") +
  geom_smooth(method = "lm", se = FALSE, aes(color = dataset), linetype = "solid", show.legend = FALSE, linewidth = 1.2) +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))
agelm
ggsave("Figures_Tables/plot_lm_age-vs-prob.png", plot = agelm, width = 6, height = 4)

# Greenness vs Classification Probability
greenlm <- ggplot(combined_data, aes(x = greenIndex, y = prob_predicted, color = dataset)) +
  geom_point(size = 0.9) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = FALSE, aes(color = dataset), linetype = "solid", show.legend = FALSE, linewidth = 1.2) +
  labs(title = NULL,
       x = "Green Index",
       y = "Classification Probability") +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))
greenlm
ggsave("Figures_Tables/lm_green_combined.png", plot = greenlm, width = 6, height = 4)

# Green Index vs  Age with Polynomial and Linear Regression
age_gi <- ggplot(combined_data, aes(x = Age, y = greenIndex, color = dataset)) +  # Correct color mapping
  geom_point(size = 0.9) +  # Points
  geom_smooth(method = "lm", se = FALSE, aes(color = dataset), linetype = "solid") +  # Linear regression line
  scale_color_manual(values = colors) +  # Dataset colors
  geom_smooth(method = "lm", se = FALSE, aes(color = dataset), linetype = "solid", show.legend = FALSE, linewidth = 1.2) +  # Linear regression overlay
  labs(title = NULL,
       x = "Age (years)",
       y = "Green Index") +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))
age_gi
# Save the plot
ggsave("Figures_Tables/lm_age_vs_greenIndex_poly_vs_linear.png", plot = age_gi, width = 6, height = 4)

# Nearest Taxon Distance vs Classification Probability
ntdlm <- ggplot(combined_data, aes(x = ntd, y = prob_predicted, color = dataset)) +
  geom_point(size = 0.9) +
  geom_smooth(method = "lm", se = FALSE, aes(color = dataset), linetype = "solid") +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = FALSE, aes(color = dataset), linetype = "solid", show.legend = FALSE, linewidth = 1.2) +
  labs(title = NULL,
       x = "Nearest Taxon Distance (M years)",
       y = "Classification Probability") +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))
ntdlm
ggsave("Figures_Tables/lm_ntd_combined.png", plot = ntdlm, width = 6, height = 4)



#Combine plots with annotations
fig7 <- (agelm + greenlm + age_gi + ntdlm) + 
  plot_layout(ncol = 2, nrow = 2, guides = "collect") +  # Arrange plots in a grid
  plot_annotation(tag_levels = "A") &  # Automatically tag plots with A, B, C, D
  theme(
    plot.tag = element_text(size = 9),  # Customize tag appearance
    legend.position = "bottom"
  )

# Save the combined plot to a single PNG file
ggsave("Figures_Tables/Fig8-regressions.png", plot = fig7, width = 8.5, height = 7)




#'------------------------------------------------------------------------------
#' @Logistic_regressions-Table4_TableS3
#-------------------------------------------------------------------------------

logmodel <- glm(correct ~ hasGlue + leafDamage + doyOfCollection + specimenQuality + tissueDevelopmentalStage + greenIndex + Age + ntd + leafKg_m2, 
             data = classProbs, 
             family = binomial)
summary(logmodel)

# Save summary to a text file
sink("result_logregression_summary.txt")
print(summary(logmodel))
sink()

# Save coefficients to a CSV file
coefficients_df <- as.data.frame(summary(logmodel)$coefficients)
write.csv(coefficients_df, "Table4_logregression_coefficients.csv", row.names = TRUE)

########### Random Forest, variable importance
library(randomForest)

classProbs$correct <- as.factor(classProbs$correct)

model_rf <- randomForest(correct ~ hasGlue + leafDamage + specimenQuality + tissueDevelopmentalStage + greenIndex + Age + ntd + leafKg_m2, data = classProbs, 
                         importance = TRUE)

# Check variable importance
importance(model_rf)

# Extract importance data from the random forest model
importance_data <- as.data.frame(importance(model_rf))

# Add a column for variable names
importance_data$Variable <- rownames(importance_data)

# Reorder columns to have Variable names first
importance_data <- importance_data[, c("Variable", "MeanDecreaseAccuracy", "MeanDecreaseGini")]

# Write the data frame to a CSV file
write.csv(importance_data, file = "TableS3_variable_importance_rf.csv", row.names = FALSE)


#'------------------------------------------------------------------------------
#' @Other-plots-and-regressions-SI
#-------------------------------------------------------------------------------

## Fig S5: Confusion matrix from coefficient-based predictions.

# Generate CM from predictions
pred_cm <- confusion_matrices_prediction_plsda(predictions_table_ntd)

# Plot
# Extract and clean confusion matrix
confusion_matrix <- as.matrix(pred_cm$mean_confusion_matrix)
rownames(confusion_matrix) <- confusion_matrix[, 1]  # Set row names
confusion_matrix <- confusion_matrix[, -1]           # Remove "Reference" column
confusion_matrix <- apply(confusion_matrix, 2, as.numeric)  # Convert to numeric

# Calculate mean accuracy
mean_accuracy <- sum(diag(pred_cm$confusion_matrix)) / sum(pred_cm$confusion_matrix)  # from performance table

# Melt the confusion matrix to long format
conf_mean_long <- tidyr::pivot_longer(
  data = pred_cm$mean_confusion_matrix,
  cols = -Reference,
  names_to = "Prediction",
  values_to = "Mean"
)

# Round Mean values to the nearest whole number
conf_mean_long$Mean <- round(conf_mean_long$Mean)

# Remove 0 values
conf_mean_long <- conf_mean_long[conf_mean_long$Mean > 0, ]

# Ensure species names match the tree tip labels
simplify_genus <- function(factor_levels) {
  unname(sapply(factor_levels, function(level) {
    parts <- strsplit(level, "_")[[1]]  # Split at the underscore
    if (length(parts) == 2) {
      paste0(substr(parts[1], 1, 3), ". ", parts[2])  # Add a period and a space, then concatenate
    } else {
      level  # Return the original level if there's no underscore
    }
  }))
}

# Apply simplify_genus to transform the values of Reference and Prediction
conf_mean_long$Reference <- simplify_genus(as.character(conf_mean_long$Reference))
conf_mean_long$Prediction <- simplify_genus(as.character(conf_mean_long$Prediction))

# Define the color palette from white to orange
cols <- colorRampPalette(c('white', '#fe9929'))

# Create the confusion matrix heatmap
cmp_plot <- ggplot(data = conf_mean_long, aes(x = Prediction, y = Reference, fill = Mean)) +
  geom_tile(color = "white") +  # White border for tiles
  geom_text(aes(label = Mean), color = '#542788', size = 4) +  # Purple text for rounded Mean values
  theme_minimal() +  # Minimal theme
  scale_fill_gradientn(colours = cols(10), limits = c(0, max(conf_mean_long$Mean, na.rm = TRUE))) +  # Gradient from white to orange
  guides(fill = "none") +  # Remove the legend
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic", size = 10),
    axis.text.y = element_text(face = "italic", size = 10), 
    panel.grid.major = element_line(color = "gray90"),  # Light major grid lines
    panel.grid.minor = element_line(color = "gray95"),  # Light minor grid lines
    panel.border = element_blank()  # Remove border
  ) +
  ggtitle(paste("Prediction from Coefficients - Mean Accuracy:", round(mean_accuracy, 3))) +  # Dynamically add mean accuracy to title
  labs(y = "True Identity", x = "Predicted Identity")

# Save the combined plot
ggsave(paste0(output_path, "confusion_matrix_prediction_sp25.png", sep = ""), plot = cmp_plot, width = 6.5, height = 6.5, dpi=300)


# Fig S7: Phylogenetic distance of predicted class against probability of classification for misclassifications.
incorrect_classifications$label <- "Incorrect Predictions"
# Create the plot
logPDproblm <- ggplot(incorrect_classifications, aes(x = log(predicted_class_PD), y = prob_predicted)) +
  geom_point(aes(color = label), size = 2, show.legend = TRUE) +  # Add legend for points
  geom_smooth(method = "lm", se = FALSE, color = "#FDE725FF") +  # Regression line matching dot color
  labs(
    title = "Classification probability (misclassifications)\ndecays with phylogenetic distance",
    x = "Phylogenetic Distance to Predicted Taxon (log(M years))",
    y = "Classification Probability"
  ) +
  scale_color_manual(values = c("Incorrect Predictions" = "#FDE725FF")) +  # Custom color for legend
  theme_minimal(base_family = "Arial") +
  theme(
    legend.title = element_blank(),  # No legend title
    axis.title = element_text(size = 10, color = "black"),  # Adjust axis title size and color
    axis.text = element_text(size = 8, color = "black"),  # Adjust axis text size and color
    plot.title = element_text(size = 12, face = "bold", color = "black", hjust = 0.5),  # Title styling
    panel.background = element_rect(fill = "gray70", color = NA),  # Black background
    plot.background = element_rect(fill = "gray70", color = NA),  # Black plot area
    panel.grid.major = element_line(color = "black"),  # Grid lines
    panel.grid.minor = element_line(color = "gray40"),
    legend.position = "bottom" # Minor grid lines
  )
logPDproblm
ggsave("Figures_Tables/Fig-S_lm_PDlog-vs-Prob_incorrect.png", plot = logPDproblm, width = 6, height = 4, dpi = 300)


# Fig. S8: Regression of LMA against classification probabilities
LMAlm <- ggplot(combined_data, aes(x = leafKg_m2, y = prob_predicted, color = dataset)) +
  geom_point(size = 0.9) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", show.legend = FALSE, linewidth = 1.2) +
  labs(title = NULL,
       x = "Leaf mass per area (Kg m-2)",
       y = "Classification Probability") +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))
LMAlm
ggsave("Figures_Tables/lm_leafLMA.png", plot = LMAlm, width = 6, height = 4)


# Fig. S9: Linear and polynomial regressions of collection Julian day against classification probabilities.
poly_doylm <- ggplot(combined_data, aes(x =  doyOfCollection, y = prob_predicted, color = dataset)) +
  geom_point(size = 0.9) +  # Points
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linetype = "dashed", size = 1.2, show.legend = FALSE) +  # Polynomial regression line
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 1.2) +  # Linear regression line
  scale_color_manual(values = colors) +  # Dataset colors
  labs(title = "Day Collected vs Classification Probability",
       subtitle = "Linear vs Polynomial Regression (Degree 2)",
       x = "Julian Calendar Day (DOY)",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
poly_doylm

# Save the plot
ggsave("Figures_Tables/lm_DOY_combined_poly_vs_linear.png", plot = poly_doylm, width = 6, height = 4)



#'------------------------------------------------------------------------------
#' @Scripts-NOT-USED
#-------------------------------------------------------------------------------

# Log Nearest Taxon Distance vs Classification Probability
logntdlm <- ggplot(combined_data, aes(x = log(ntd), y = prob_predicted, color = dataset)) +
  geom_point(size = 0.9) +
  geom_smooth(method = "lm", se = FALSE, aes(color = dataset), linetype = "solid") +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = FALSE, aes(color = dataset), linetype = "solid", show.legend = FALSE, linewidth = 1.2) +
  labs(title = "log(NTD) vs Classification Probability",
       x = "Nearest Taxon Distance (log(M years))",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
logntdlm
ggsave("herbarium-predictors-analysis/lm_LogNTD_combined.png", plot = logntdlm, width = 6, height = 4)



### Proportion-correct-by-NTD

#predictions_table_ntd_pd <- fread(paste0(output_path, "predictions_table_plsda25spp.csv"))

predictions <- predictions_table_ntd_pd

# Calculate proportions for each unique NTD
ntd_summary <- predictions[, .(
  Total = .N,
  Correct = sum(correct == TRUE),
  Incorrect = sum(correct == FALSE)
), by = ntd]

# Add proportions
ntd_summary[, `:=`(
  Correct_Proportion = Correct / Total,
  Incorrect_Proportion = Incorrect / Total
)]

# Reshape data for plotting
plot_data <- melt(
  ntd_summary,
  id.vars = "ntd",
  measure.vars = c("Correct_Proportion", "Incorrect_Proportion"),
  variable.name = "Classification",
  value.name = "Proportion"
)

# Update classification labels
plot_data[, Classification := ifelse(Classification == "Correct_Proportion", "Correct", "Incorrect")]

# Filter plot_data to include only correct classifications
correct_data <- plot_data[Classification == "Correct"]

# Generate scatterplot for correct classifications only
ntd_scatter_plot <- ggplot(correct_data, aes(x = ntd, y = Proportion)) +
  geom_point(size = 1.5, color = "#21908CFF") +  # Scatter dots with size 0.8 and blue color
  geom_smooth(method = "lm", se = FALSE, color = "#21908CFF", linewidth = 1) +  # Linear regression line
  labs(
    title = "Proportion of Correct Classifications by NTD",
    x = "Nearest Taxon Distance (NTD)",
    y = "Proportion"
  ) +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"  # Remove legend as only one classification is shown
  )

# Display the plot
print(ntd_scatter_plot)

# Save the scatterplot to a PNG file
ggsave(paste0(output_path, "plot_Pcorrect_by_ntd.png"), plot = ntd_scatter_plot, width = 6.5, height = 4)

#'------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------
