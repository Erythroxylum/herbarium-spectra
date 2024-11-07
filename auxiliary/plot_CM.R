#'------------------------------------------------------------------------------
#' @title Plotting Confusion Matrices
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

#-------------------------------------------------------------------------------
#' @Plot-Confusion-Matrices
#-------------------------------------------------------------------------------

# Load data
CM_data <- CM_plsda_testing
# set data name for plot header
file_name <- deparse(substitute(CM_plsda_testing))

# Extract mean matrix Melt the confusion matrix to long format
conf_mean_long <- melt(CM_data$mean_confusion_matrix, id.vars = "Reference", variable.name = "Prediction", value.name = "Mean")

# Round Mean values to the nearest whole number
conf_mean_long$Mean <- round(conf_mean_long$Mean)

# Remove 0 values
conf_mean_long <- conf_mean_long[conf_mean_long$Mean > 0, ] 

# Define the color palette from white to orange
cols <- colorRampPalette(c('white', '#fe9929'))

# Plot using ggplot with Mean as the fill
cm_plot <- ggplot(data = conf_mean_long, aes(x = Prediction, y = Reference, fill = Mean)) +
  geom_tile(color = "white") +  # White border for tiles
  geom_text(aes(label = Mean), color = '#542788', size = 4) +  # Purple text for rounded Mean values
  theme_minimal() +  # Minimal theme
  scale_fill_gradientn(colours = cols(10), limits = c(0, max(conf_mean_long$Mean, na.rm = TRUE))) +  # Gradient from white to orange
  guides(fill = "none") +  # Remove the legend
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic", size = 12),
    axis.text.y = element_text(face = "italic", size = 12),
    panel.grid.major = element_line(color = "gray90"),  # Light major grid lines
    panel.grid.minor = element_line(color = "gray95"),  # Light minor grid lines
    panel.border = element_blank()  # Remove border
  ) +
  ggtitle(paste(file_name)) +  # Dynamically add file name to title
  labs(y = "True Identity", x = "Predicted Identity")

cm_plot

# Save the plot
ggsave(paste(file_name, ".pdf", sep=""), plot = cm_plot, width = 7.5, height = 7.5)



