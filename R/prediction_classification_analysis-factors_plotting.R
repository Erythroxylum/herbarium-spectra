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

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results
root_path <- getwd()
output_path <- paste0(root_path, "/Figures_Tables/")

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Read data
frame <- fread(paste0(root_path, 
                      "/data/dataHUH2024_sp25leaf561_ref5nm_450-2400.csv"))

#-------------------------------------------------------------------------------
#' @Data_reshape  
#-------------------------------------------------------------------------------

# Get files from meta data, traits, and spectra.

meta <- frame[, c("collector", "accession", "accession_leaf", "leaf", "scan","species",
                  "ddmmyyScanned", "doy", "absoluteAge", "herbQuality",
                  "damage", "glue", "leafStage", "greenIndex")]

traits <- frame[, c("leafKg_m2", "leafThickness")]

meta$sample <- 1:nrow(meta)

species <- frame$species
species <- sub(" ", "_", species)

spectra <- frame[, .SD, .SDcols = 23:ncol(frame)]


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
phylo <- read.tree("herbarium-predictors-analysis/phylogram_pd_TimeTree5.tre")

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
correct_classifications <- predictions %>%
  filter(true_class == predicted_class) %>%
  mutate(dataset = "Correct Predictions")

# Generate data frame with incorrect classifications
incorrect_classifications <- predictions %>%
  filter(true_class != predicted_class) %>%
  mutate(dataset = "Incorrect Predictions")

# Add a label for the full dataset
classProbs <- predictions %>%
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

# Plot for herbQualityity
p1 <- ggplot(cor_inc_data, aes(x = herbQuality, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  labs(x = "Specimen Quality", y = "Predicted Probabilities") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF", "#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Plot for glue
p2 <- ggplot(cor_inc_data, aes(x = glue, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  labs(x = "Glue", y = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF", "#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Plot for damage
p3 <- ggplot(cor_inc_data, aes(x = damage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  labs(x = "Damage Status", y = "Predicted Probabilities") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF", "#FDE725FF")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Plot for leafStage
p4 <- ggplot(cor_inc_data, aes(x = leafStage, y = prob_predicted, fill = dataset)) +
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
# herbQuality
# Filter data for Correct Predictions
correct_data <- cor_inc_data[dataset == "Correct Predictions"]

# Plot for Specimen Quality within Correct Predictions
p1c <- ggplot(correct_data, aes(x = herbQuality, y = prob_predicted, fill = dataset)) +
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
p1i <- ggplot(incorrect_data, aes(x = herbQuality, y = prob_predicted, fill = dataset)) +
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
p2_correct <- ggplot(correct_glue, aes(x = glue, y = prob_predicted, fill = dataset)) +
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
p2_incorrect <- ggplot(incorrect_glue, aes(x = glue, y = prob_predicted, fill = dataset)) +
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
# Damage

# Filter data for Correct Predictions
correct_damage <- cor_inc_data[dataset == "Correct Predictions"]

# Filter data for Incorrect Predictions
incorrect_damage <- cor_inc_data[dataset == "Incorrect Predictions"]

# Plot for Damage (Correct Predictions)
p3_correct <- ggplot(correct_damage, aes(x = damage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("major", "medium"), c("major","minor"),c("major","none"),c("medium","minor"),c("medium","none"),c("minor", "none")),  # Replace with actual levels of damage
    map_signif_level = TRUE,
    test = "t.test"
  ) +
  labs(x = "Damage Status", y = "Predicted Probabilities", title = "Damage: Correct Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, by = 0.02))

p3_correct <- ggplot(correct_damage, aes(x = damage, y = prob_predicted, fill = dataset)) +
  geom_boxplot(outlier.colour = "black") +
  geom_signif(
    comparisons = list(c("major", "medium"), c("major", "minor"), c("major", "none"),
                       c("medium", "minor"), c("medium", "none"), c("minor", "none")),  # Levels of damage
    map_signif_level = TRUE,
    test = "t.test",
    y_position = c(0.122, 0.126, 0.132, 0.112, 0.136, 0.122)  # Adjust y-positions to prevent overlap
  ) +
  labs(x = "Damage Status", y = "Predicted Probabilities", title = "Damage: Correct Predictions") +
  theme_minimal() +
  scale_fill_manual(values = c("#21908CFF")) +
  scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, by = 0.02))


# Plot for Damage (Incorrect Predictions)
p3_incorrect <- ggplot(incorrect_damage, aes(x = damage, y = prob_predicted, fill = dataset)) +
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


###### 1. Test for herbQualityity
kruskal_test_herbQualityity <- kruskal.test(prob_predicted ~ herbQuality, data = classProbs)
print(kruskal_test_herbQualityity)

# Pairwise comparisons using Dunn's test
dunn_herbQuality <- dunn.test(classProbs$prob_predicted, classProbs$herbQuality, method = "bh")
dunn_herbQuality
# [1] "good - medium" "good - poor"   "medium - poor"
# p 4.535799e-01 4.527629e-06 4.197685e-06

######  2. Test for glue
kruskal_test_glue <- kruskal.test(prob_predicted ~ glue, data = correct_classifications)
print(kruskal_test_glue)
# Kruskal-Wallis chi-squared = 3.7555, df = 1, p-value = 0.05263

#######  3. Test for damage
kruskal_test_damage <- kruskal.test(prob_predicted ~ damage, data = correct_classifications)
print(kruskal_test_damage)
#Kruskal-Wallis chi-squared = 15.65, df = 3, p-value = 0.001338

# Pairwise comparisons
dunn_damage <- dunn.test(correct_classifications$prob_predicted, correct_classifications$damage, method = "bh")
print(dunn_damage$P.adjusted)
#$P.adjusted
#[1] 0.1505601144 0.2545721306 0.1375234946 0.1282708638 0.2344388191 0.0005490144
#$comparisons
#[1] "major - medium" "major - minor"  "medium - minor" "major - none"   "medium - none" 
#[6] "minor - none"  

#######  4. Test for leafStage
kruskal_test_leafStage <- kruskal.test(prob_predicted ~ leafStage, data = correct_classifications)
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
agelm <- ggplot(combined_data, aes(x = absoluteAge, y = prob_predicted, color = dataset)) +
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
ggsave("herbarium-predictors-analysis/plot_lm_age-vs-prob.png", plot = agelm, width = 6, height = 4)

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
ggsave("herbarium-predictors-analysis/lm_green_combined.png", plot = greenlm, width = 6, height = 4)

# Green Index vs  Age with Polynomial and Linear Regression
age_gi <- ggplot(combined_data, aes(x = absoluteAge, y = greenIndex, color = dataset)) +  # Correct color mapping
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
ggsave("herbarium-predictors-analysis/lm_age_vs_greenIndex_poly_vs_linear.png", plot = age_gi, width = 6, height = 4)

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
ggsave("herbarium-predictors-analysis/lm_ntd_combined.png", plot = ntdlm, width = 6, height = 4)



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

logmodel <- glm(correct ~ glue + damage + doy + herbQuality + leafStage + greenIndex + absoluteAge + ntd + leafKg_m2, 
             data = classProbs, 
             family = binomial)
summary(logmodel)

# Save summary to a text file
sink("herbarium-predictors-analysis/result_logregression_summary.txt")
print(summary(logmodel))
sink()

# Save coefficients to a CSV file
coefficients_df <- as.data.frame(summary(logmodel)$coefficients)
write.csv(coefficients_df, "herbarium-predictors-analysis/Table4_logregression_coefficients.csv", row.names = TRUE)

########### Random Forest, variable importance
library(randomForest)

classProbs$correct <- as.factor(classProbs$correct)

model_rf <- randomForest(correct ~ glue + damage + herbQuality + leafStage + greenIndex + absoluteAge + ntd + leafKg_m2, data = classProbs, 
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
write.csv(importance_data, file = "herbarium-predictors-analysis/TableS3_variable_importance_rf.csv", row.names = FALSE)


#'------------------------------------------------------------------------------
#' @Other-plots-and-regressions-SI
#-------------------------------------------------------------------------------

## Fig S5: Confusion matrix from coefficient-based predictions.

# Generate CM from predictions
pred_cm <- confusion_matrices_prediction_plsda(predictions_table)

# Plot
# Extract and clean confusion matrix
confusion_matrix <- as.matrix(pred_cm$mean_confusion_matrix)
rownames(confusion_matrix) <- confusion_matrix[, 1]  # Set row names
confusion_matrix <- confusion_matrix[, -1]           # Remove "Reference" column
confusion_matrix <- apply(confusion_matrix, 2, as.numeric)  # Convert to numeric

# Calculate mean accuracy
mean_accuracy <- sum(diag(pred_cm$confusion_matrix)) / sum(pred_cm$confusion_matrix)  # from performance table

# Melt the confusion matrix to long format
conf_mean_long <- melt(pred_cm$mean_confusion_matrix, id.vars = "Reference", variable.name = "Prediction", value.name = "Mean")

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
ggsave("herbarium-predictors-analysis/lm_leafLMA.png", plot = LMAlm, width = 6, height = 4)


# Fig. S9: Linear and polynomial regressions of collection Julian day against classification probabilities.
poly_doylm <- ggplot(combined_data, aes(x = doy, y = prob_predicted, color = dataset)) +
  geom_point(size = 0.9) +  # Points
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linetype = "dashed", size = 1.2, show.legend = FALSE) +  # Polynomial regression line
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", size = 1.2) +  # Linear regression line
  scale_color_manual(values = colors) +  # Dataset colors
  labs(title = "Day Collected vs Classification Probability",
       subtitle = "Linear vs Polynomial Regression (Degree 2)",
       x = "Julian Calendar Day (DOY)",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
poly_doylm

# Save the plot
ggsave("herbarium-predictors-analysis/lm_DOY_combined_poly_vs_linear.png", plot = poly_doylm, width = 6, height = 4)



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



### PLOT probability by scan number 0,1,2

#predictions_table_ntd_pd <- fread(paste0(output_path, "predictions_table_plsda25spp.csv"))

# Restrict scans to the first three per accession_leaf and renumber rows
predictions_scan <- predictions_table_ntd_pd[, .SD[1:3], by = accession_leaf]
predictions_scan[, scan_index := seq_len(.N) - 1, by = accession_leaf]  # Renumber rows as 0, 1, 2

# Ensure that only groups with all three scans are included
valid_predictions <- predictions_scan[, if (.N == 3) .SD, by = accession_leaf]

# Extract probability of the predicted class
valid_predictions[, predicted_probability := mapply(function(class, row) {
  row[[class]]  # Directly access the column matching the predicted class
}, predicted_class, split(.SD, seq_len(nrow(.SD))), SIMPLIFY = TRUE),
.SDcols = colnames(predictions_scan)[grepl("^[A-Z]", colnames(predictions_scan))]]

# Prepare data for plotting (include accession_leaf for validation)
plot_data <- valid_predictions[, .(accession_leaf, scan_index, correct, predicted_probability)]

# Remove accession_leaf groups where all correct values are the same
filtered_plot_data <- plot_data[, if (!(all(correct) || all(!correct))) .SD, by = accession_leaf]

## Chi-squared test
# Create contingency table: Count correct and incorrect predictions for each scan
contingency_table <- filtered_plot_data[, .N, by = .(scan_index, correct)]
contingency_table <- dcast(contingency_table, scan_index ~ correct, value.var = "N", fill = 0)
setnames(contingency_table, c("scan_index", "FALSE", "TRUE"), c("scan_index", "Incorrect", "Correct"))

# Perform Chi-Square Test
chi_square_test <- chisq.test(contingency_table[, .(Correct, Incorrect)])

# Print results
print("Contingency Table:")
print(contingency_table)
print("Chi-Square Test Results:")
print(chi_square_test)

# Add a column for total predictions per scan
contingency_table[, Total := Correct + Incorrect]

# Calculate proportions
contingency_table[, Correct_Proportion := Correct / Total]
contingency_table[, Incorrect_Proportion := Incorrect / Total]

# Melt the data for plotting proportions
proportion_data <- melt(
  contingency_table,
  id.vars = "scan_index",
  measure.vars = c("Incorrect_Proportion", "Correct_Proportion"),
  variable.name = "Classification",
  value.name = "Proportion"
)

# Update classification labels for the legend
proportion_data[, Classification := ifelse(Classification == "Correct_Proportion", "Correct", "Incorrect")]

# Create labels for incorrect and correct counts
contingency_table[, Label := paste0("Incorrect: ", Incorrect, "\nCorrect: ", Correct)]

# Reorder levels of correct to ensure "Correct" is on the left
plot_data[, correct := factor(correct, levels = c(TRUE, FALSE))]

# Generate the probabilities boxplot
boxplot <- ggplot(plot_data, aes(x = factor(scan_index), y = predicted_probability, fill = correct)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7) +
  labs(
    title = "Prediction Probabilities",
    x = "Scan",
    y = "Probability of Predicted Class",
    fill = "Classification"
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#21908CFF", "FALSE" = "#FDE725FF"),
    labels = c("TRUE" = "Correct", "FALSE" = "Incorrect")  # Custom legend labels
  ) +
  #scale_y_continuous(limits = c(0.035, 0.125)) +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  )
# Save the plot
#ggsave("classification_probabilities_boxplot_with_counts.png", plot = boxplot, width = 8, height = 6)

# Combine the plots side by side using patchwork
combined_plot <- bar_chart + boxplot + plot_layout(ncol = 2)

# Display the combined plot
print(combined_plot)

# Save the combined plot to a single PDF
ggsave(paste0(output_path, "plot_pred_by_scan-Pvalue0.2119.png"), plot = combined_plot, width = 8.5, height = 4)


#'------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------