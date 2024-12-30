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

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("auxiliary/pred_coef.R")

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
coef <- fread(paste0(root_path, "/out_classify/coefficients_plsda.csv"))
coef$Predictor <- gsub("`", "", coef$Predictor)

# Get class names from coefficients
classes <- colnames(coef)[!colnames(coef) %in% c("Predictor", "Iteration")]

# Predict probabilities
predictions <- predict_plsda(
  spectra = spectra,
  coefficients = coef,
  classes = classes,
  true_classes = species,  # True class for validation
  meta = meta,  # Metadata table
  output_file = "plsda_probabilities_with_metadata.csv"
)


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
predictions[, True_Class := gsub(" ", "_", True_Class)]
predictions[, Predicted_Class := gsub(" ", "_", Predicted_Class)]

# Calculate NTD for each True_Class
ntd <- distance_matrix_dt[species1 != species2, .(
  NTD = min(PD_value, na.rm = TRUE)
), by = species1]

# Rename `species1` to `True_Class` for merging compatibility
setnames(ntd, "species1", "True_Class")

# Merge NTD into predictions
predictions <- merge(predictions, ntd, by = "True_Class", all.x = TRUE)

# Add the phylogenetic distance to the Predicted_Class
predictions <- merge(
  predictions, 
  distance_matrix_dt,
  by.x = c("True_Class", "Predicted_Class"),
  by.y = c("species1", "species2"),
  all.x = TRUE
)

# Rename the distance column for clarity
setnames(predictions, "PD_value", "Predicted_Class_PD")

# Fill `Predicted_Class_PD` for correct classifications with 0
predictions[Correct_Prediction == TRUE, Predicted_Class_PD := 0]

##save
saveRDS(predictions, "herbarium-predictors-analysis/huh_ref5nm_25spp_classification.rds")


#'------------------------------------------------------------------------------
#' @Plot_Prob_by_scan
#-------------------------------------------------------------------------------

predictions <- readRDS("herbarium-predictors-analysis/huh_ref5nm_25spp_classification.rds")

# Ensure the data.table structure is correct
predictions_scan <- as.data.table(predictions)

# Restrict scans to the first three per accession_leaf and renumber rows
predictions_scan <- predictions_scan[, .SD[1:3], by = accession_leaf]
predictions_scan[, Scan_Index := seq_len(.N) - 1, by = accession_leaf]  # Renumber rows as 0, 1, 2

# Ensure that only groups with all three scans are included
valid_predictions <- predictions_scan[, if (.N == 3) .SD, by = accession_leaf]

# Extract probability of the predicted class
valid_predictions[, Predicted_Probability := mapply(function(class, row) {
  row[[paste0("prob_", class)]]
}, Predicted_Class, split(.SD, seq_len(nrow(.SD))), SIMPLIFY = TRUE),
.SDcols = patterns("^prob_")]

# Prepare data for plotting (include accession_leaf for validation)
plot_data <- valid_predictions[, .(accession_leaf, Scan_Index, Correct_Prediction, Predicted_Probability)]

# Remove accession_leaf groups where all Correct_Prediction values are the same
filtered_plot_data <- plot_data[, if (!(all(Correct_Prediction) || all(!Correct_Prediction))) .SD, by = accession_leaf]

## Chi-squared test
# Create contingency table: Count correct and incorrect predictions for each scan
contingency_table <- filtered_plot_data[, .N, by = .(Scan_Index, Correct_Prediction)]
contingency_table <- dcast(contingency_table, Scan_Index ~ Correct_Prediction, value.var = "N", fill = 0)
setnames(contingency_table, c("Scan_Index", "FALSE", "TRUE"), c("Scan_Index", "Incorrect", "Correct"))

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
  id.vars = "Scan_Index",
  measure.vars = c("Incorrect_Proportion", "Correct_Proportion"),
  variable.name = "Classification",
  value.name = "Proportion"
)

# Update classification labels for the legend
proportion_data[, Classification := ifelse(Classification == "Correct_Proportion", "Correct", "Incorrect")]

# Create labels for incorrect and correct counts
contingency_table[, Label := paste0("Incorrect: ", Incorrect, "\nCorrect: ", Correct)]

# Generate the bar plot
bar_chart <- ggplot(proportion_data, aes(x = factor(Scan_Index), y = Proportion, fill = Classification)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(
    data = contingency_table,
    aes(x = factor(Scan_Index), y = 0.69, label = Label),  # Adjust y to place labels
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    title = "Proportion of Correct and Incorrect Predictions by Scan",
    x = "Scan",
    y = "Proportion",
    fill = "Classification"
  ) +
  scale_fill_manual(values = c("Incorrect" = "red", "Correct" = "blue")) +
  scale_y_continuous(limits = c(0, 0.72)) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  )

# Reorder levels of Correct_Prediction to ensure "Correct" is on the left
plot_data[, Correct_Prediction := factor(Correct_Prediction, levels = c(TRUE, FALSE))]

# Generate the probabilities boxplot
boxplot <- ggplot(plot_data, aes(x = factor(Scan_Index), y = Predicted_Probability, fill = Correct_Prediction)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7) +
  labs(
    title = "Prediction Probabilities",
    x = "Scan",
    y = "Probability of Predicted Class",
    fill = "Classification"
  ) +
  scale_fill_manual(
    values = c("TRUE" = "blue", "FALSE" = "red"),
    labels = c("TRUE" = "Correct", "FALSE" = "Incorrect")  # Custom legend labels
  ) +
  #scale_y_continuous(limits = c(0.035, 0.125)) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  )
# Save the plot
#ggsave("classification_probabilities_boxplot_with_counts.pdf", plot = boxplot, width = 8, height = 6)

# Combine the plots side by side using patchwork
combined_plot <- bar_chart + boxplot + plot_layout(ncol = 2)

# Save the combined plot to a single PDF
ggsave("classification_by_scan_combined_plots.pdf", plot = combined_plot, width = 8.5, height = 4)

# Display the combined plot
print(combined_plot)



#'------------------------------------------------------------------------------
#' @Proportion-correct-by-NTD
#-------------------------------------------------------------------------------

### Seems to be a problem where several taxa were not predicted very well?

predictions <- readRDS("herbarium-predictors-analysis/huh_ref5nm_25spp_classification.rds")

# Calculate proportions for each unique NTD
ntd_summary <- predictions[, .(
  Total = .N,
  Correct = sum(Correct_Prediction == TRUE),
  Incorrect = sum(Correct_Prediction == FALSE)
), by = NTD]

# Add proportions
ntd_summary[, `:=`(
  Correct_Proportion = Correct / Total,
  Incorrect_Proportion = Incorrect / Total
)]

# Reshape data for plotting
plot_data <- melt(
  ntd_summary,
  id.vars = "NTD",
  measure.vars = c("Correct_Proportion", "Incorrect_Proportion"),
  variable.name = "Classification",
  value.name = "Proportion"
)

# Update classification labels
plot_data[, Classification := ifelse(Classification == "Correct_Proportion", "Correct", "Incorrect")]

# Filter plot_data to include only correct classifications
correct_data <- plot_data[Classification == "Correct"]

# Generate the line plot for correct classifications only
line_plot <- ggplot(correct_data, aes(x = NTD, y = Proportion, color = Classification, group = Classification)) +
  geom_line(linewidth = 1, color = "blue") +  # Set color to blue for "Correct"
  geom_point(size = 2, color = "blue") +      # Set point color to blue
  labs(
    title = "Proportion of Correct Classifications by NTD",
    x = "Nearest Taxon Distance (NTD)",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"  # Remove legend as only one classification is shown
  )

# Save the plot
ggsave("proportion_by_ntd_no_bins.pdf", plot = line_plot, width = 8, height = 6)

# Display the plot
print(line_plot)


#'------------------------------------------------------------------------------
#' @Subset_DataFrame
#-------------------------------------------------------------------------------

# Generate data frame with correct classifications
correct_classifications <- predictions %>%
  filter(True_Class == Predicted_Class) %>%
  mutate(dataset = "Correct Predictions")

# Generate data frame with incorrect classifications
incorrect_classifications <- predictions %>%
  filter(True_Class != Predicted_Class) %>%
  mutate(dataset = "Incorrect Predictions")

# Add a label for the full dataset
classProbs <- predictions %>%
  mutate(dataset = "All")


#'------------------------------------------------------------------------------
#' @Plots-Prob_vs_predictors
#-------------------------------------------------------------------------------

#Correct classification bar plots
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr) # for ggarrange function
library(dunn.test)
library(ggsignif)

# Pick datasets
all_data <- bind_rows(classProbs, correct_classifications, incorrect_classifications)
all_data <- bind_rows(correct_classifications, incorrect_classifications)

# Set y-axis limits for consistency
ymin = 0
ymax = 0.2

# Plot for herbQuality
p1 <- ggplot(all_data, aes(x = herbQuality, y = Prob, fill = dataset)) +
  geom_boxplot(outlier.colour = "red") +
  labs(x = "Herb Quality", y = "Predicted Probabilities") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen", "lightcoral")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Plot for glue
p2 <- ggplot(all_data, aes(x = glue, y = Prob, fill = dataset)) +
  geom_boxplot(outlier.colour = "red") +
  labs(x = "Glue Type", y = "Predicted Probabilities") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen", "lightcoral")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Plot for damage
p3 <- ggplot(all_data, aes(x = damage, y = Prob, fill = dataset)) +
  geom_boxplot(outlier.colour = "red") +
  labs(x = "Damage Status", y = "Predicted Probabilities") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen", "lightcoral")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Plot for leafStage
p4 <- ggplot(all_data, aes(x = leafStage, y = Prob, fill = dataset)) +
  geom_boxplot(outlier.colour = "red") +
  labs(x = "Leaf Stage", y = "Predicted Probabilities") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen", "lightcoral")) +
  scale_y_continuous(limits = c(ymin, ymax))

# Arrange all plots in a grid
plots4 <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# Save the plot as a PDF
pdf("../herbarium_spectra_results/leaf_plsdav1/PLSDA-regressions_correct-incorrect.pdf", width = 10, height = 10)
print(plots4)
dev.off()


###### 1. Test for herbQuality
kruskal_test_herbQuality <- kruskal.test(Prob ~ herbQual, data = classProbs)
print(kruskal_test_herbQuality)

# Pairwise comparisons using Dunn's test
dunn_herbQual <- dunn.test(classProbs$Prob, classProbs$herbQual, method = "bh")
dunn_herbQual
# [1] "good - medium" "good - poor"   "medium - poor"
# p 4.535799e-01 4.527629e-06 4.197685e-06

######  2. Test for glue
kruskal_test_glue <- kruskal.test(Prob ~ glue, data = correct_classifications)
print(kruskal_test_glue)
# Kruskal-Wallis chi-squared = 3.7555, df = 1, p-value = 0.05263

#######  3. Test for damage
kruskal_test_damage <- kruskal.test(Prob ~ damage, data = correct_classifications)
print(kruskal_test_damage)
#Kruskal-Wallis chi-squared = 15.65, df = 3, p-value = 0.001338

# Pairwise comparisons
dunn_damage <- dunn.test(correct_classifications$Prob, correct_classifications$damage, method = "bh")
print(dunn_damage$P.adjusted)
#$P.adjusted
#[1] 0.1505601144 0.2545721306 0.1375234946 0.1282708638 0.2344388191 0.0005490144
#$comparisons
#[1] "major - medium" "major - minor"  "medium - minor" "major - none"   "medium - none" 
#[6] "minor - none"  

#######  4. Test for leafStage
kruskal_test_leafStage <- kruskal.test(Prob ~ leafStage, data = correct_classifications)
print(kruskal_test_leafStage)
# Kruskal-Wallis chi-squared = 0.0793, df = 1, p-value = 0.78

dev.off()


#'------------------------------------------------------------------------------
#' @Linear_Regressions_vs_numeric_predictors
#-------------------------------------------------------------------------------

# Combine datasets for plotting
combined_data <- bind_rows(classProbs, correct_classifications, incorrect_classifications)

# Absolute Age vs Classification Probability
agelm <- ggplot(combined_data, aes(x = absoluteAge, y = Prob, color = dataset)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Absolute Age vs Classification Probability",
       x = "Absolute Age",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
agelm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_age_combined.pdf", plot = agelm, width = 6, height = 4)

# Greenness vs Classification Probability
greenlm <- ggplot(combined_data, aes(x = greenIndex, y = Prob, color = dataset)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Greenness vs Classification Probability",
       x = "Green Index",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
greenlm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_green_combined.pdf", plot = greenlm, width = 6, height = 4)

# Nearest Taxon Distance vs Classification Probability
logntdlm <- ggplot(combined_data, aes(x = log(ntd), y = Prob, color = dataset)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "NTD vs Classification Probability",
       x = "Nearest Taxon Distance (M years)",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
logntdlm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_LogNTD_combined.pdf", plot = logntdlm, width = 6, height = 4)

# Nearest Taxon Distance vs Classification Probability
ntdlm <- ggplot(combined_data, aes(x = ntd, y = Prob, color = dataset)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "NTD vs Classification Probability",
       x = "Nearest Taxon Distance (M years)",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
ntdlm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_ntd_combined.pdf", plot = ntdlm, width = 6, height = 4)

# Log Nearest Taxon Distance vs Classification Probability
logntdlm <- ggplot(combined_data, aes(x = log(ntd), y = Prob, color = dataset)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "log(NTD) vs Classification Probability",
       x = "Nearest Taxon Distance (log(M years)",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
logntdlm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_LogNTD_combined.pdf", plot = logntdlm, width = 6, height = 4)

# Phylogenetic Distance vs Classification Probability
PDproblm <- ggplot(incorrect_classifications, aes(x = PD_value, y = Prob)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Predicted PD (incorrect classificaitons) vs Classification Probability",
       x = "Phylogenetic distance to predicted taxon (M years)",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
PDproblm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_PD-vs-Prob_incorrect.pdf", plot = PDproblm, width = 4, height = 4)

# Log Phylogenetic Distance vs Classification Probability
logPDproblm <- ggplot(incorrect_classifications, aes(x = log(PD_value), y = Prob)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Log Predicted PD (incorrect classificaitons)\n vs Classification Probability",
       x = "Phylogenetic distance to predicted taxon (log(M years)",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
logPDproblm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_PDlog-vs-Prob_incorrect.pdf", plot = logPDproblm, width = 4, height = 4)

#'------------------------------------------------------------------------------
#' @Linear_Regressions_deeper
#-------------------------------------------------------------------------------

# Regress classification probability on absoluteAge (numeric)
lm_age <- lm(Prob ~ absoluteAge, data = correct_classifications)
summary(lm_age)
plot(lm_age)

#'------------------------------------------------------------------------------
#' @Logistic_regressions
#-------------------------------------------------------------------------------

model <- glm(correct ~ glue + damage + herbQual + leafStage + greenIndex + absoluteAge + ntd, 
             data = classProbs, 
             family = binomial)
summary(model)


########### Random Forest, variable importance
library(randomForest)

model_rf <- randomForest(correct ~ glue + damage + herbQual + leafStage + greenIndex + absoluteAge + ntd, data = classProbs, importance = TRUE)

# Check variable importance
importance(model_rf)

# Extract importance data from the random forest model
importance_data <- as.data.frame(importance(model_rf))

# Add a column for variable names
importance_data$Variable <- rownames(importance_data)

# Reorder columns to have Variable names first
importance_data <- importance_data[, c("Variable", "MeanDecreaseAccuracy", "MeanDecreaseGini")]

# Write the data frame to a CSV file
write.csv(importance_data, file = "../herbarium_spectra_results/leaf_plsdav1/variable_importance_rf.csv", row.names = FALSE)


#'------------------------------------------------------------------------------
#' @Age-vs-GreenIndex
#-------------------------------------------------------------------------------

plot(frame$absoluteAge, frame$greenIndex,
     main = "Linear Regression of Green Index vs. Absolute Age",
     xlab = "Absolute Age",
     ylab = "Green Index",
     pch = 16,              # Set point shape
     col = "blue")          # Set point color

# Add a linear regression line
abline(lm(greenIndex ~ absoluteAge, data = frame), col = "red", lwd = 2)

# Optional: Add regression line equation and R-squared to the plot
model <- lm(greenIndex ~ absoluteAge, data = frame)
eq <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                 list(a = format(coef(model)[1], digits = 2), 
                      b = format(coef(model)[2], digits = 2), 
                      r2 = format(summary(model)$r.squared, digits = 3)))
mtext(as.expression(eq), 3, line = -1.5, col="red")








#######################################################
################## old
## Extract the probabilities of incorrect predictions.
#This line compares the predicted classes (plsClasses_all) to the true species labels (herb_spec_all$meta$SpShort).
#outputs logical vector where TRUE indicates a misclassification (the prediction is different from the true species) and FALSE indicates a correct classification.
#incorrect_classifications <- plsClasses_all != herb_spec$meta$SpShort
incorrect_classifications <- predicted_classes != species

# Get the probability of the incorrect class for those misclassifications
incorrect_probs <- plsProbs_all[[1]][cbind(1:nrow(plsProbs_all[[1]]), match(predicted_classes, colnames(plsProbs_all[[1]])))]

# Convert to a data frame and add the column name
incorrect_probs <- data.frame(IncorrectProb = incorrect_probs)

# Set correct classifications to NA
incorrect_probs[incorrect_classifications == FALSE] <- NA  

# Combine incorrect probabilities with metadata
meta_probs <- cbind(meta, IncorrectProb = incorrect_probs)

# Convert absoluteAge to numeric, if needed
meta_probs$absoluteAge <- as.numeric(as.character(meta_probs$absoluteAge))






#####################################################################
#####################################################################
## Summing incorrect class probabilities

##########################################################
# Summing incorrect classification probabilities
incorrect_sum <- sum(incorrect_probs, na.rm = TRUE)
print(paste("Total incorrect classification probability:", incorrect_sum))

##########################################################
# Linear regression model
herb_data$damage <- herb_spec$meta$damage
lm_model <- lm(Prob ~ absoluteAge + herbQual + glue + damage + leafStage, data = correct_classifications)
summary(lm_model)

# Logistic regression model (converting probabilities into binary outcome for logistic regression)
# Define threshold (e.g., misclassification probabilities above 0.5 are considered "incorrect")
threshold <- 0.5
herb_data$IncorrectClass <- ifelse(herb_data$IncorrectProb > threshold, 1, 0)

# Fit logistic regression
log_model <- glm(IncorrectClass ~ absoluteAge + herbQuality + glue + damage + leafStage, 
                 data = herb_data, family = "binomial")
summary(log_model)



##########################################################
# Visualization - Logistic model

# Plot predicted probabilities from the logistic model
herb_data$predicted_logit <- predict(log_model, newdata = herb_data, type = "response")

ggplot(herb_data, aes(x = absoluteAge, y = predicted_logit)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "purple") +
  labs(title = "Logistic Regression: Predicted Probability of Incorrect Classification",
       x = "Absolute Age",
       y = "Predicted Probability of Incorrect Classification")

# Optionally, plot logistic regression by herbQuality
ggplot(herb_data, aes(x = herbQuality, y = predicted_logit)) +
  geom_boxplot() +
  labs(title = "Logistic Regression: Herbarium Quality vs Predicted Incorrect Classification Probability",
       x = "Herbarium Quality",
       y = "Predicted Probability")




############ ANOVA
anova_glue <- aov(Prob ~ glue, data = correct_classifications)
summary(anova_glue)

