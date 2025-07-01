#'------------------------------------------------------------------------------
#' @title Plotting functions for leaf models
#'------------------------------------------------------------------------------

#' @description A script to plot the output from 'Plotting functions for leaf
#' leaf models'
#'
#' @return several plots

#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(ape)
library(patchwork)
library(phytools)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------


#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

getwd()

#-------------------------------------------------------------------------------
#' @Classification-Performance-Table-5
#-------------------------------------------------------------------------------

# Define the file paths and metadata

files <- list(
  list(path = "out_classify/performance_lda_testing_sp10.rds", dataset = "Herbarium", rank = "species", model = "LDA", samples = "10spp"),
  list(path = "out_classify/performance_plsda_testing_sp10.rds", dataset = "Herbarium", rank = "species", model = "PLS-DA", samples = "10spp"),
  list(path = "out_classify_kot/performance_lda_testing_sp10.rds", dataset = "Pressed", rank = "species", model = "LDA", samples = "10spp"),
  list(path = "out_classify_kot/performance_plsda_testing_sp10.rds", dataset = "Pressed", rank = "species", model = "PLS-DA", samples = "10spp"),
  list(path = "out_classify_genus/performance_lda_testing_sp10.rds", dataset = "Herbarium", rank = "genus", model = "LDA", samples = "6genera"),
  list(path = "out_classify_genus/performance_plsda_testing_sp10.rds", dataset = "Herbarium", rank = "genus", model = "PLS-DA", samples = "6genera"),
  list(path = "out_classify_kot_LatinGenus/performance_lda_testing_sp10.rds", dataset = "Pressed", rank = "genus", model = "LDA", samples = "6genera"),
  list(path = "out_classify_kot_LatinGenus/performance_plsda_testing_sp10.rds", dataset = "Pressed", rank = "genus", model = "PLS-DA", samples = "6genera"),
  list(path = "out_classify/performance_lda_testing.rds", dataset = "Herbarium", rank = "species", model = "LDA", samples = "25spp"),
  list(path = "out_classify/performance_plsda_testing.rds", dataset = "Herbarium", rank = "species", model = "PLS-DA", samples = "25spp"),
  list(path = "out_classify_genus/performance_lda_testing.rds", dataset = "Herbarium", rank = "genus", model = "LDA", samples = "17genera"),
  list(path = "out_classify_genus/performance_plsda_testing.rds", dataset = "Herbarium", rank = "genus", model = "PLS-DA", samples = "17genera")
   )

# Initialize an empty data frame for combined statistics
combined_stats <- data.frame()

# Loop through the files and process each
for (file_info in files) {
  # Load the performance data
  perf <- readRDS(file_info$path)
  
  # Convert to data frame
  perf_df <- as.data.frame(perf$performance)
  
  # Check if ncomp exists in the data frame
  has_ncomp <- "ncomp" %in% colnames(perf_df)
  
  # Calculate summary statistics
  summary_stats <- perf_df %>%
    summarise(
      mean_accuracy = mean(Accuracy, na.rm = TRUE),
      sd_accuracy = sd(Accuracy, na.rm = TRUE),
      mean_precision = mean(Precision, na.rm = TRUE),
      sd_precision = sd(Precision, na.rm = TRUE),
      mean_balanced_accuracy = mean(`Balanced Accuracy`, na.rm = TRUE),
      sd_balanced_accuracy = sd(`Balanced Accuracy`, na.rm = TRUE),
      Ncomponents = if (has_ncomp) mean(ncomp, na.rm = TRUE) else NA_real_  # Handle missing ncomp
    ) %>%
    mutate(
      Dataset = file_info$dataset,
      Rank = file_info$rank,
      Model = file_info$model,
      Samples = file_info$samples
    )
  
  # Append to the combined stats table
  combined_stats <- bind_rows(combined_stats, summary_stats)
}

# Select and arrange columns in the specified order
final_table <- combined_stats %>%
  select(Rank, Dataset, Model, Samples, mean_accuracy, sd_accuracy, mean_precision, sd_precision, mean_balanced_accuracy, sd_balanced_accuracy, Ncomponents)

# Combine mean and SD into a single column
final_table <- combined_stats %>%
  mutate(
    Accuracy = sprintf("%.2f ± %.2f", mean_accuracy, sd_accuracy),
    Precision = sprintf("%.2f ± %.2f", mean_precision, sd_precision),
    `Balanced Accuracy` = sprintf("%.2f ± %.2f", mean_balanced_accuracy, sd_balanced_accuracy)
  ) %>%
  select(Dataset, Rank, Model, Samples,Ncomponents, Accuracy, Precision, `Balanced Accuracy`, )


# Save the final table to a CSV file
write.csv(final_table, "Figures_Tables/Table3_classification_performance_maxAcc.csv", row.names = FALSE)



#-------------------------------------------------------------------------------
#' @Plot-Confusion-Matrices_PLSDA-25sp_Fig5
#-------------------------------------------------------------------------------

## load phylogeny
# Load classification confusion matrix output
CM_data <- readRDS("out_classify_kot/CM_plsda_testing_sp10.rds")

# Read in Tree
phylo <- read.tree("herbarium-predictors-analysis/phylogram_pd_TimeTree5.tre")

# Set data name for plot header
file_name <- "Pressed-leaf, PLS-DA Validation"

# Extract and clean confusion matrix
confusion_matrix <- as.matrix(CM_data$mean_confusion_matrix)
rownames(confusion_matrix) <- confusion_matrix[, 1]  # Set row names
confusion_matrix <- confusion_matrix[, -1]           # Remove "Reference" column
confusion_matrix <- apply(confusion_matrix, 2, as.numeric)  # Convert to numeric

# Calculate mean accuracy
# Assuming CM_data is loaded
# Extract the mean_confusion_matrix
confusion_matrix <- CM_data$mean_confusion_matrix

# Convert the confusion matrix to a numeric matrix (excluding the 'Reference' column)
confusion_matrix_numeric <- as.matrix(confusion_matrix[, -1])
# Extract the diagonal values (correct classifications)
correct_classifications <- diag(confusion_matrix_numeric)
# Calculate the mean accuracy
mean_accuracy <- mean(correct_classifications)
# Print the mean accuracy
print(paste("Mean Accuracy:", round(mean_accuracy, 2)))

#manually
#mean_accuracy <- 0.743  # from performance table

# Melt the confusion matrix to long format
conf_mean_long <- melt(CM_data$mean_confusion_matrix, id.vars = "Reference", variable.name = "Prediction", value.name = "Mean")

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

# Extract the actual tip order from the ladderized tree
tree_plot <- ggtree(phylo)
tree_tip_order <- tree_plot$data %>%
  dplyr::filter(isTip) %>%  # Ensure dplyr is loaded
  dplyr::arrange(y) %>%
  dplyr::pull(label)

tree_tip_order <- simplify_genus(as.character(tree_tip_order))

# Reorder the Reference and Prediction factors based on the tree tip order
conf_mean_long$Reference <- factor(conf_mean_long$Reference, levels = tree_tip_order)
conf_mean_long$Prediction <- factor(conf_mean_long$Prediction, levels = tree_tip_order)

# Define the color palette from white to orange
cols <- colorRampPalette(c('white', '#fe9929'))

# Create the confusion matrix heatmap
cm_plot <- ggplot(data = conf_mean_long, aes(x = Prediction, y = Reference, fill = Mean)) +
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
  #ggtitle(paste(file_name, "- Mean Accuracy:", round(mean_accuracy, 3))) +  # Dynamically add mean accuracy to title
  labs(y = NULL, x = "Predicted Identity")



################################
### Plot phylogram 

# Define time scale (assume branch lengths are in millions of years)
tree_plot$data$x <- max(tree_plot$data$x, na.rm = TRUE) - tree_plot$data$x
max_time <- max(tree_plot$data$x, na.rm = TRUE)
#time_breaks <- seq(0, max_time, by = 50)
time_breaks <- seq(max_time, 0, by = -50)
time_labels <- c(0, 50, 100, 150, 200, 250, 300, 350, 400)

# Add the phylogeny plot with branch lengths and x-axis for years
tree_plot <- ggtree(phylo) +
  geom_vline(xintercept = time_breaks, color = "gray80", linetype = "dashed") +  # Add vertical bars
  scale_x_continuous(breaks = time_breaks, labels = time_labels) +  # Ascending time scale
  theme_tree2() +  # Use theme_tree2 to enable axis customization
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(size=10, vjust = 20),
    axis.text.y = element_blank(),  # Remove y-axis text for the tree
    axis.title.y = element_text(size = 10, angle = 90, hjust = 0.5),  # Add y-axis label
    panel.grid = element_blank()
  ) +
  labs(y = "True Identity", x = "M Years")  # Add axis labels


######
##Combine the phylogeny and confusion matrix using patchwork
combined_plot <- tree_plot + cm_plot + plot_layout(widths = c(1.5, 4)) 

# Save the combined plot
ggsave(paste("Figures_Tables/Fig6_", file_name, " with phylogram.png", sep = ""), plot = combined_plot, width = 10, height = 7, dpi=300)



#-------------------------------------------------------------------------------
#' @PLSDA_VIP_Plots_by_species-Fig-S3
#-------------------------------------------------------------------------------

vip_plsda <- fread("out_classify/varImp_plsda.csv")

# Convert the data table to a data frame
df <- as.data.frame(vip_plsda)

# Convert to data.table
df <- as.data.table(df)

# Step 1: Rename "Variable" to "Wavelength"
setnames(df, "Variable", "Wavelength")

# Step 2: Remove backticks from "Wavelength" and convert to numeric
df[, Wavelength := as.numeric(gsub("`", "", Wavelength))]

# Step 3: Melt the data for long format
df_long <- melt(df, id.vars = c("Model", "Wavelength"), 
                variable.name = "Species", 
                value.name = "VIP_Value")

# Modify factor levels for Reference and Prediction
simplify_genus <- function(factor_levels) {
  sapply(factor_levels, function(level) {
    parts <- strsplit(level, "_")[[1]]  # Split at the underscore
    if (length(parts) == 2) {
      paste0(substr(parts[1], 1, 3), ". ", parts[2])  # Add a period and a space, then concatenate
    } else {
      level  # Return the original level if there's no underscore
    }
  })
}

# Apply simplify_genus to transform the values of Reference and Prediction
df_long$Species <- factor(simplify_genus(as.character(df_long$Species )))

# Step 4: Calculate mean, 50% and 90% quantiles, and CV for each Species at each Wavelength
vip_stats <- df_long[, .(
  Mean_VIP = mean(VIP_Value, na.rm = TRUE),
  low_50 = quantile(VIP_Value, 0.25, na.rm = TRUE),
  high_50 = quantile(VIP_Value, 0.75, na.rm = TRUE),
  low_90 = quantile(VIP_Value, 0.05, na.rm = TRUE),
  high_90 = quantile(VIP_Value, 0.95, na.rm = TRUE),
  CV_VIP = ifelse(mean(VIP_Value, na.rm = TRUE) == 0, NA, sd(VIP_Value, na.rm = TRUE) / mean(VIP_Value, na.rm = TRUE)) # not informative
), by = .(Wavelength, Species)]

vip_stats_save <- vip_stats

# Plotting
vipsp <- ggplot(vip_stats, aes(x = Wavelength)) +
  # 90% interval ribbon
  geom_ribbon(aes(ymin = low_90, ymax = high_90, fill = "90% Range"), alpha = 0.5) +
  # 50% interval ribbon
  geom_ribbon(aes(ymin = low_50, ymax = high_50, fill = "50% Range"), alpha = 0.9) +
  # Mean VIP line
  geom_line(aes(y = Mean_VIP, color = "Mean VIP"), linewidth = 0.3) +
  labs(
    title = "Variable Importance in Projection: Herbarium Leaf PLSDA Model",
    x = "Wavelength (nm)",
    y = "VIP Values",
    fill = "Intervals",  # Legend title for ribbons
    color = "Metric"     # Legend title for lines
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.5),
    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size = 10)
  ) +
  # Scale x-axis to have ticks every 200 nm
  scale_x_continuous(breaks = seq(0, max(vip_stats$Wavelength, na.rm = TRUE), by = 200)) +
  scale_y_continuous() +
  # Custom colors for ribbons and lines
  scale_fill_manual(values = c("90% Range" = "blue", "50% Range" = "blue")) +
  scale_color_manual(values = c("Mean VIP" = "black")) +
  # Facet for each species using the italicized labels, arranged in 4 columns
  facet_wrap(~ Species, scales = "free_y", ncol = 4)

# Save the plot
ggsave("Figures_Tables/Fig-S3_VIP-herb-plsda-25spp.png", plot = vipsp, width = 8.5, height = 9.5, dpi=300)



#-------------------------------------------------------------------------------
#' @PLSDA_VIP_Plots_across_datasets-Fig-S4
#-------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(gridExtra)

# Define the file paths and metadata
file_info <- list(
  list(path = "out_classify/varImp_plsda_sp10.csv", dataset = "Herbarium", rank = "species", model = "PLS-DA", samples = "10 spp"),
  list(path = "out_classify_kot/varImp_plsda_sp10.csv", dataset = "Pressed", rank = "species", model = "PLS-DA", samples = "10 spp"),
  list(path = "out_classify_genus/varImp_plsda_sp10.csv", dataset = "Herbarium", rank = "genus", model = "PLS-DA", samples = "6 genera"),
  list(path = "out_classify_kot_LatinGenus/varImp_plsda_sp10.csv", dataset = "Pressed", rank = "genus", model = "PLS-DA", samples = "6 genera"),
  list(path = "out_classify/varImp_plsda.csv", dataset = "Herbarium", rank = "species", model = "PLS-DA", samples = "25 spp"),
  list(path = "out_classify_genus/varImp_plsda.csv", dataset = "Herbarium", rank = "genus", model = "PLS-DA", samples = "17 genera")
)

# Function to create a plot for a given file
create_vip_plot <- function(file_info) {
  # Load data
  vip_data <- fread(file_info$path)
  
  # Rename "Variable" to "Wavelength" and remove backticks
  setnames(vip_data, "Variable", "Wavelength")
  vip_data[, Wavelength := as.numeric(gsub("`", "", Wavelength))]
  
  # Melt the data for long format
  vip_long <- melt(vip_data, id.vars = c("Model", "Wavelength"), 
                   variable.name = "Species", 
                   value.name = "VIP_Value")
  
  # Calculate statistics
  vip_stats <- vip_long[, .(
    Mean_VIP = mean(VIP_Value, na.rm = TRUE),
    low_50 = quantile(VIP_Value, 0.25, na.rm = TRUE),
    high_50 = quantile(VIP_Value, 0.75, na.rm = TRUE),
    low_90 = quantile(VIP_Value, 0.05, na.rm = TRUE),
    high_90 = quantile(VIP_Value, 0.95, na.rm = TRUE),
    CV_VIP = ifelse(mean(VIP_Value, na.rm = TRUE) == 0, NA, 
                    sd(VIP_Value, na.rm = TRUE) / mean(VIP_Value, na.rm = TRUE))
  ), by = .(Wavelength)]
  
  # Create the plot
  p <- ggplot(vip_stats, aes(x = Wavelength)) +
    geom_ribbon(aes(ymin = low_90, ymax = high_90), fill = "blue", alpha = 0.5) +
    geom_ribbon(aes(ymin = low_50, ymax = high_50), fill = "blue", alpha = 0.9) +
    geom_line(aes(y = Mean_VIP), color = "black", linewidth = 0.8) +
    #geom_line(aes(y = CV_VIP, color = "red"), linewidth = 0.8) + # uniformly too low to be informative
    labs(title = paste(file_info$dataset, file_info$model, file_info$samples),
         y = "VIP Values", x = NULL) +  # Remove x-axis label
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray", linewidth = 0.5),
      panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    scale_x_continuous(breaks = seq(0, max(vip_stats$Wavelength, na.rm = TRUE), by = 200)) +
    # Custom colors for ribbons and lines
    scale_fill_manual(values = c("90% Range" = "lightblue", "50% Range" = "blue")) +
    scale_color_manual(values = c("Mean VIP" = "black"))
  
  return(p)
}



# Generate all plots
plots <- lapply(file_info, create_vip_plot)

# Adjust x-axis labels for the bottom plots on each page
for (i in c(6)) { # c(4, 6)
  plots[[i]] <- plots[[i]] + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                                   axis.ticks.x = element_line(),
                                   axis.title.x = element_text(size = 10)) +
    labs(x = "Wavelength (nm)")
}

# Arrange and save the plots, 7 x 9 layout
pdf("Figures_Tables/Fig-S4_VIP-pressed-herbv2.pdf", width = 7, height = 9)
# First page with six plots
grid.arrange(grobs = plots[1:6], ncol = 1)
dev.off()

##

###############################
#### Find peaks

# Function to identify peaks in VIP data
find_peaks <- function(vip_stats) {
  # Ensure data is sorted by Wavelength
  vip_stats <- vip_stats[order(Wavelength)]
  
  # Identify peaks: A point is a peak if its Mean_VIP is greater than its neighbors
  peaks <- vip_stats[
    (shift(Mean_VIP, type = "lag", fill = -Inf) < Mean_VIP) & 
      (shift(Mean_VIP, type = "lead", fill = -Inf) < Mean_VIP)
  ]
  
  return(peaks[, .(Wavelength, Mean_VIP)])
}

# Example: Find peaks for one specific dataset
# Assuming vip_stats for one dataset is available after plotting
file_to_analyze <- file_info[[1]]  # Replace with the desired file_info index
vip_data <- fread(file_to_analyze$path)

# Rename "Variable" to "Wavelength" and remove backticks
setnames(vip_data, "Variable", "Wavelength")
vip_data[, Wavelength := as.numeric(gsub("`", "", Wavelength))]

# Melt the data for long format
vip_long <- melt(vip_data, id.vars = c("Model", "Wavelength"), 
                 variable.name = "Species", 
                 value.name = "VIP_Value")

# Calculate statistics
vip_stats <- vip_long[, .(
  Mean_VIP = mean(VIP_Value, na.rm = TRUE),
  low_50 = quantile(VIP_Value, 0.25, na.rm = TRUE),
  high_50 = quantile(VIP_Value, 0.75, na.rm = TRUE),
  low_90 = quantile(VIP_Value, 0.05, na.rm = TRUE),
  high_90 = quantile(VIP_Value, 0.95, na.rm = TRUE)
), by = .(Wavelength)]

# Find peaks in the Mean_VIP values
peaks <- find_peaks(vip_stats)

# Print the peaks
print(peaks)

# Optionally save the peaks to a CSV file
fwrite(peaks, "Figures_Tables/Fig-S4_VIP_Peaks_Specific_Dataset.csv")

#'------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------

