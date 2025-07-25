#-------------------------------------------------------------------------------
# Function to generate and summarize confusion matrices (mean and SD)

# Load necessary libraries
library(data.table)
library(dplyr)
library(tidyr)
library(caret)  # For confusionMatrix function

#-------------------------------------------------------------------------------
# PLSDA

confusion_matrices_plsda <- function(meta_split, species_split, spectra_split, models, ncomp, threads = 1) {
  
  # Data for predictions
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"
  
  # Ensure `species_split` is a factor
  species_split <- factor(species_split)
  
  # Initialize list to store confusion matrices
  confusion_matrices <- vector("list", length(models))
  
  # Loop over each model to generate confusion matrices
  for (x in 1:length(models)) {
    # Predict with the current model
    predictions <- predict(models[[x]], 
                           newdata = frame, 
                           ncomp = ncomp)
    
    # Ensure predictions and species_split have the same levels
    predictions <- factor(predictions, 
                          levels = levels(species_split))
    
    # Generate the confusion matrix
    herb_conmat <- confusionMatrix(predictions, species_split)
    
    # Convert the confusion matrix table to a data frame
    pcm_d <- as.data.frame(herb_conmat$table)
    
    # Filter out any rows with NA values in Prediction or Reference
    pcm_d <- pcm_d %>% filter(!is.na(Prediction) & !is.na(Reference))
    
    # Calculate percentages
    pcm_d$RefSum <- unlist(lapply(pcm_d$Reference, function(x) sum(pcm_d$Freq[pcm_d$Reference == x])))
    pcm_d$RefPer <- round(pcm_d$Freq / pcm_d$RefSum * 100, digits = 0)
    
    # Store the processed confusion matrix in the list
    confusion_matrices[[x]] <- pcm_d
  }
  
  # Combine all confusion matrices for summarization
  combined_results <- do.call(rbind, confusion_matrices)
  
  # Calculate mean and standard deviation frequencies for each prediction/reference pair
  summary_matrix <- combined_results %>%
    group_by(Prediction, Reference) %>%
    summarise(
      MeanPer = mean(RefPer, na.rm = TRUE),
      SDPer = sd(RefPer, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Pivot to create wide-format confusion matrices for means and standard deviations
  mean_confusion_matrix <- summary_matrix %>%
    select(Prediction, Reference, MeanPer) %>%
    pivot_wider(names_from = Prediction, values_from = MeanPer, values_fill = 0)
  
  sd_confusion_matrix <- summary_matrix %>%
    select(Prediction, Reference, SDPer) %>%
    pivot_wider(names_from = Prediction, values_from = SDPer, values_fill = 0)
  
  # Return both matrices
  return(list(mean_confusion_matrix = mean_confusion_matrix, sd_confusion_matrix = sd_confusion_matrix))
}

#-------------------------------------------------------------------------------
# CM from predictions (pred_coef_plsda.R)

confusion_matrices_prediction_plsda <- function(output_dt) {
  library(caret)
  library(dplyr)
  library(tidyr)
  
  # Ensure true and predicted are factors with same levels
  all_levels <- union(unique(output_dt$true_class), unique(output_dt$predicted_class))
  true_classes <- factor(output_dt$true_class, levels = all_levels)
  predicted_classes <- factor(output_dt$predicted_class, levels = all_levels)
  
  # Generate the confusion matrix
  cm <- confusionMatrix(predicted_classes, true_classes)
  
  # Extract the confusion matrix table as a data frame
  cm_table <- as.data.frame(cm$table)
  colnames(cm_table) <- c("Reference", "Prediction", "Frequency")
  
  # Calculate percentages for each reference class
  cm_table <- cm_table %>%
    group_by(Reference) %>%
    mutate(Percentage = round((Frequency / sum(Frequency)) * 100, 1))
  
  # Pivot to wide format
  mean_confusion_matrix <- cm_table %>%
    select(Reference, Prediction, Percentage) %>%
    pivot_wider(names_from = Prediction, values_from = Percentage, values_fill = 0)
  
  return(list(
    confusion_matrix = cm$table,
    cm_table = cm_table,
    mean_confusion_matrix = mean_confusion_matrix
  ))
}

#-------------------------------------------------------------------------------
# LDA

confusion_matrices_lda <- function(meta_split, species_split, spectra_split, models, threads = 1) {
  
  # Data for predictions
  frame <- cbind(species_split, spectra_split)
  colnames(frame)[1] <- "species"
  
  # Ensure `species_split` is a factor
  species_split <- factor(species_split)
  
  # Initialize list to store confusion matrices
  confusion_matrices <- vector("list", length(models))
  
  # Loop over each model to generate confusion matrices
  for (x in 1:length(models)) {
    # Predict with the current model
    predictions <- predict(models[[x]], newdata = frame)
    
    # Ensure predictions and species_split have the same levels
    predictions <- factor(predictions, levels = levels(species_split))
    
    # Generate the confusion matrix
    herb_conmat <- confusionMatrix(predictions, species_split)
    
    # Convert the confusion matrix table to a data frame
    pcm_d <- as.data.frame(herb_conmat$table)
    
    # Filter out any rows with NA values in Prediction or Reference
    pcm_d <- pcm_d %>% filter(!is.na(Prediction) & !is.na(Reference))
    
    # Calculate percentages
    pcm_d$RefSum <- unlist(lapply(pcm_d$Reference, function(x) sum(pcm_d$Freq[pcm_d$Reference == x])))
    pcm_d$RefPer <- round(pcm_d$Freq / pcm_d$RefSum * 100, digits = 0)
    
    # Store the processed confusion matrix in the list
    confusion_matrices[[x]] <- pcm_d
  }
  
  # Combine all confusion matrices for summarization
  combined_results <- do.call(rbind, confusion_matrices)
  
  # Calculate mean and standard deviation frequencies for each prediction/reference pair
  summary_matrix <- combined_results %>%
    group_by(Prediction, Reference) %>%
    summarise(
      MeanPer = mean(RefPer, na.rm = TRUE),
      SDPer = sd(RefPer, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Pivot to create wide-format confusion matrices for means and standard deviations
  mean_confusion_matrix <- summary_matrix %>%
    select(Prediction, Reference, MeanPer) %>%
    pivot_wider(names_from = Prediction, values_from = MeanPer, values_fill = 0)
  
  sd_confusion_matrix <- summary_matrix %>%
    select(Prediction, Reference, SDPer) %>%
    pivot_wider(names_from = Prediction, values_from = SDPer, values_fill = 0)
  
  # Return both matrices
  return(list(mean_confusion_matrix = mean_confusion_matrix, sd_confusion_matrix = sd_confusion_matrix))
}

