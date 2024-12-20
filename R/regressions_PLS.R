#'------------------------------------------------------------------------------
#' @title Regressions of classification probabilities against herbarium quality 
#' predictors
#'------------------------------------------------------------------------------
#' 
#' @description Script to assemble classification probabilities, phylogenetic 
#' distances and herbarium quality metadata and analyze under regressions 
#' and ANCOVA
#' 
#' @return plots
#' 
#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------

library(spectrolab)
library(caret)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr) # for ggarrange function
library(dunn.test)
library(ggsignif)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------
setwd("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/results/LDA/")

# Select the root_folder to read data and export results
root_path <- getwd()

#'------------------------------------------------------------------------------
#' @Load-PLS-Predictions 
#-------------------------------------------------------------------------------

##########################################################
## Load PLS model predictions

testing <- fread(paste0(root_path,"/results/HUH/ref/HUH_450-2400/pls_HUH_testing_obs-pred.csv"))
training <- fread(paste0(root_path,"/results/HUH/ref/HUH_450-2400/pls_HUH_training_obs-pred.csv"))
predict_data <- rbind(testing,training)


#'------------------------------------------------------------------------------
#' @Plots-Prob_vs_predictors
#-------------------------------------------------------------------------------

#Correct classification bar plots
# Load necessary libraries

# Set y-axis limits for consistency
ymin = 0
ymax = 0.2

### herbQual, good to poor P 0.051
boxplot(predict_data$predicted_sd ~ predict_data$herbQuality)
boxplot(predict_data$predicted_se ~ predict_data$herbQuality) # same
boxplot(abs(predict_data$residuals_mean) ~ predict_data$herbQuality)
boxplot(predict_data$residuals_sd ~ predict_data$herbQuality)
boxplot(predict_data$residuals_se ~ predict_data$herbQuality) # same

### Damage nope
boxplot(predict_data$predicted_sd ~ predict_data$damage)
boxplot(abs(predict_data$residuals_mean) ~ predict_data$damage)
boxplot(predict_data$residuals_sd ~ predict_data$damage)

### Glue, very significant for sd, but in wrong direction. Sig for residuals P 0.01 !!!
boxplot(predict_data$predicted_sd ~ predict_data$glue)
boxplot(abs(predict_data$residuals_mean) ~ predict_data$glue)
boxplot(predict_data$residuals_sd ~ predict_data$glue)

### leaf Stage, SD sig 0.02, residuals sig but wrong way
boxplot(predict_data$predicted_sd ~ predict_data$leafStage)
boxplot(abs(predict_data$residuals_mean) ~ predict_data$leafStage)
boxplot(predict_data$residuals_sd ~ predict_data$leafStage)


## predicted SD
kruskal_test_herbQuality <- kruskal.test(predict_data$predicted_sd ~ predict_data$leafStage, data = predict_data) # not sig
print(kruskal_test_herbQuality)
dunn_herbQual <- dunn.test(predict_data$predicted_sd, predict_data$leafStage, method = "bh")
dunn_herbQual # good to poor P 0.051

## Residuals
kruskal_test_herbQuality <- kruskal.test(abs(predict_data$residuals_mean )~ predict_data$leafStage, data = predict_data) # not sig
print(kruskal_test_herbQuality)

dunn_herbQual <- dunn.test(abs(predict_data$residuals_mean), predict_data$herbQual, method = "bh")
dunn_herbQual # not sig

## Residuals SD
kruskal_test_herbQuality <- kruskal.test(predict_data$predicted_se ~ predict_data$leafStage, data = predict_data) # not sig
print(kruskal_test_herbQuality)
dunn_herbQual <- dunn.test(predict_data$predicted_se, predict_data$herbQuality, method = "bh")
dunn_herbQual # not sig



### PLOTS
# Plot for herbQuality
p1 <- ggplot(predict_data, aes(x = predict_data$herbQual, y = predict_data$predicted_se, fill = predict_data$herbQual)) +
  geom_boxplot(outlier.colour = "red") +
  labs(x = "Herb Quality", y = "Predicted Probabilities") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen", "lightcoral"))
  #scale_y_continuous(limits = c(ymin, ymax))

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

# Absolute Age vs residuals, slight neg
agelm <- ggplot(predict_data, aes(x = predict_data$absoluteAge, y = abs(predict_data$residuals_mean))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Absolute Age vs Classification Probability",
       x = "Absolute Age",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
agelm

# Absolute Age vs predicted_sd, negative, reverse expectration
agelm <- ggplot(predict_data, aes(x = predict_data$absoluteAge, y = predict_data$predicted_sd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Absolute Age vs Classification Probability",
       x = "Absolute Age",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
agelm

#ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_age_combined.pdf", plot = agelm, width = 6, height = 4)


# Greenness residuals slight neg
greenlm <- ggplot(predict_data, aes(x = predict_data$greenIndex, abs(predict_data$residuals_mean))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Greenness vs Classification Probability",
       x = "Green Index",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
greenlm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_green_combined.pdf", plot = greenlm, width = 6, height = 4)

# Greenness predicted sd NEGATIVE. Green is good.
greenlm <- ggplot(predict_data, aes(x = predict_data$greenIndex, y = predict_data$predicted_sd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Greenness vs Classification Probability",
       x = "Green Index",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
greenlm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_green_combined.pdf", plot = greenlm, width = 6, height = 4)



# DOY residuals 
greenlm <- ggplot(predict_data, aes(x = predict_data$doy, y = abs(predict_data$residuals_mean))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Greenness vs Classification Probability",
       x = "Green Index",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
greenlm
#ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_green_combined.pdf", plot = greenlm, width = 6, height = 4)

# DOY predicted sd slight pos. Could be 
greenlm <- ggplot(predict_data, aes(x = predict_data$doy, y = predict_data$predicted_sd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "Greenness vs Classification Probability",
       x = "Green Index",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
greenlm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_green_combined.pdf", plot = greenlm, width = 6, height = 4)


# Nearest Taxon Distance vs Classification Probability
logntdlm <- ggplot(predict_data, aes(x = log(ntd), y = residuals_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  labs(title = "NTD vs Classification Probability",
       x = "Nearest Taxon Distance (M years)",
       y = "Classification Probability") +
  theme(legend.title = element_blank())
logntdlm
ggsave("../herbarium_spectra_results/leaf_plsdav1/lm_LogNTD_combined.pdf", plot = logntdlm, width = 6, height = 4)




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

