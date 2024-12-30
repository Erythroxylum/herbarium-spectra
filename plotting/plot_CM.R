#'------------------------------------------------------------------------------
#' @title_Plotting_Confusion_Matrices
#'------------------------------------------------------------------------------

#' @description A script to take output from leaf_plsda.R and plot matrixes
#'
#' @return .png

#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ggplot2)
library(ggtree)
library(ape)
library(data.table)
library(patchwork)
library(dplyr)

#-------------------------------------------------------------------------------
#' @Plot-Confusion-Matrices_PLSDA-25sp_Fig3
#-------------------------------------------------------------------------------

#### With phylogeny

library(ggplot2)
library(ggtree)
library(ape)
library(data.table)
library(patchwork)
library(phytools)

# Load classification confusion matrix output
CM_data <- readRDS("out_classify/CM_plsda_testing.rds")

# Read in Tree (optional)
phylo <- read.tree("herbarium-predictors-analysis/phylogram_pd_TimeTree5.tre")

# Set data name for plot header
file_name <- "Herbarium PLS-DA Validation"

# Extract and clean confusion matrix
confusion_matrix <- as.matrix(CM_data$mean_confusion_matrix)
rownames(confusion_matrix) <- confusion_matrix[, 1]  # Set row names
confusion_matrix <- confusion_matrix[, -1]           # Remove "Reference" column
confusion_matrix <- apply(confusion_matrix, 2, as.numeric)  # Convert to numeric

# Calculate mean accuracy
mean_accuracy <- 0.743  # from performance table

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
  ggtitle(paste(file_name, "- Mean Accuracy:", round(mean_accuracy, 3))) +  # Dynamically add mean accuracy to title
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
  labs(y = "True Identity", x = "Millions of Years Ago")  # Add axis labels


# Combine the phylogeny and confusion matrix using patchwork
combined_plot <- tree_plot + cm_plot + plot_layout(widths = c(1.5, 4)) 

# Save the combined plot
ggsave(paste(file_name, " with phylogram.pdf", sep = ""), plot = combined_plot, width = 10, height = 7)



################################
### Plot proportional ultrametric tree 

# (optional) Remove branch lengths for cladogram
clado <- phylo
clado$edge.length <- NULL

# Add the phylogeny plot without tip labels
tree_plot_clado <- ggtree(clado) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title.y = element_text(size = 10, angle = 90, hjust = 0.5),
    panel.grid = element_blank()
  ) +
  labs(y = "True Identity")

# Combine the phylogeny and confusion matrix using patchwork
combined_plot <- tree_plot_clado + cm_plot + plot_layout(widths = c(1.5, 4)) 

# Save the combined plot
ggsave(paste(file_name, " with cladogram.pdf", sep = ""), plot = combined_plot, width = 10, height = 7)


