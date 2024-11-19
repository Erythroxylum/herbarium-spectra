#'------------------------------------------------------------------------------
#' @title Regressions of classification probabilities against herbarium quality 
#' predictors
#'------------------------------------------------------------------------------

#' @description Script to assemble classification probabilities, phylogenetic 
#' distances and herbarium quality metadata and analyze under regressions and 
#' ANCOVA
 
#' @return 
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
#' @Working_directory
#-------------------------------------------------------------------------------
setwd("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/results/LDA/")

# Select the root_folder to read data and export results
root_path <- "C:/Users/jog4076/Downloads"
root_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/herbarium-spectra/"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Select the file of interest
frame <- fread(paste0(root_path, 
                      "../fullDataHUH2024_sp25leaf636_noResample_400-2300.csv"))

#-------------------------------------------------------------------------------
#' @Data_reshape  
#-------------------------------------------------------------------------------

# Get files from meta data, traits, and spectra.

meta <- frame[, c("collector", "accession", "accession_leaf", "leaf", "scan","species",
                  "ddmmyyScanned", "absoluteAge", "herbQuality",
                  "damage", "glue", "leafStage", "greenIndex")]

traits <- frame[, c("leafKg_m2", "leafThickness")]

meta$sample <- 1:nrow(meta)

species <- frame$species
species <- sub(" ", "_", species)

spectra <- frame[, .SD, .SDcols = 22:ncol(frame)]


#'------------------------------------------------------------------------------
#' @Build_DataFrame
#-------------------------------------------------------------------------------

##########################################################
## Load PLSDA model
models_plsda_regression <- readRDS("../leaf_plsdav1/models_plsda.rds")

## Apply the trained model to predict probabilities and class labels
plsProbs_all <- predict(models_plsda_regression[1], newdata = as.matrix(spectra), type = "prob")[[1]]  # First element in list
plsClasses_all <- predict(models_plsda_regression[1], newdata = as.matrix(spectra))[[1]]

## Retrieve true classes and metadata
true_class <- species  # Assuming 'species' is aligned with spectra rows
predicted_class <- plsClasses_all

## Combine `accession_leaf_scan` by joining columns from meta
sample <- as.character(paste(meta$accession_leaf, meta$scan, sep = "_"))

## Prepare probabilities for the true class in `Prob` column
Prob <- sapply(1:nrow(plsProbs_all), function(i) plsProbs_all[i, true_class[i]])

## Determine if each prediction is correct
correct <- true_class == predicted_class

#### Add Phylo Distance
# Load the phylogenetic diversity data
phyloDiv <- fread("../herbarium-predictors-analysis/PD-s25-scenario3.csv")

# Ensure column names match the species names in the phyloDiv data
setnames(phyloDiv, old = "V1", new = "species1")

# Reshape `phyloDiv` to long format using all columns except `species_name`
phyloDiv_long <- phyloDiv %>%
  pivot_longer(
    cols = -species1,
    names_to = "species2",
    values_to = "PD_value"
  )

## Add MNTD
# Convert data to a data table if it's not already
species_data <- as.data.table(phyloDiv) # from below

# Calculate the smallest non-zero value for each species and round to the nearest whole number
NTD <- species_data %>%
  rowwise() %>%
  mutate(ntd = round(min(c_across(-species1)[c_across(-species1) > 0], na.rm = TRUE))) %>%
  select(species1, ntd) %>%
  rename(species = species1)

# Save to a CSV file
write.csv(NTD, file = "../herbarium-predictors-analysis/MNTD.csv", row.names = FALSE)

# Rename the 'species' column in MNTD to 'true_class' for merging compatibility
NTD <- NTD %>% rename(true_class = species)

# Assemble the data frame
classProbs <- data.frame(
  sample = sample,
  true_class = as.factor(true_class),
  predicted_class = as.factor(predicted_class),
  correct = as.factor(correct),
  herbQual = as.factor(meta$herbQual),
  damage = as.factor(meta$damage),
  glue = as.factor(meta$glue),
  leafStage = as.factor(meta$leafStage),
  absoluteAge = as.numeric(meta$absoluteAge),
  greenIndex = as.numeric(meta$greenIndex),
  Prob = as.numeric(Prob)
)

# Ensure factors are correctly defined
classProbs$correct <- as.factor(ifelse(classProbs$true_class == classProbs$predicted_class, TRUE, FALSE))
classProbs$true_class <- as.factor(classProbs$true_class)
classProbs$predicted_class <- as.factor(classProbs$predicted_class)


# Merge with NTD data
classProbs <- classProbs %>%
  left_join(NTD, by = "true_class")

# Merge with PD for true_class and predicted_class
classProbs <- classProbs %>%
  left_join(phyloDiv_long, by = c("true_class" = "species1", "predicted_class" = "species2"))

# Check the structure of the updated result_table
str(classProbs)

# Save to a CSV file
saveRDS(classProbs, file = "../herbarium-predictors-analysis/FullData_classProbs.rds")

#'------------------------------------------------------------------------------
#' @Subset_DataFrame
#-------------------------------------------------------------------------------

# Generate data frame with correct classifications
correct_classifications <- classProbs %>%
  filter(true_class == predicted_class) %>%
  mutate(dataset = "Correct Predictions")

# Generate data frame with incorrect classifications
incorrect_classifications <- classProbs %>%
  filter(true_class != predicted_class) %>%
  mutate(dataset = "Incorrect Predictions")

# Add a label for the full dataset
classProbs <- classProbs %>%
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
p1 <- ggplot(all_data, aes(x = herbQual, y = Prob, fill = dataset)) +
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

