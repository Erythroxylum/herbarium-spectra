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


# Extract mean matrix and melt it into long format
#conf_sd <- CM_lda_testing$confusionMatrices_sd
conf_mean <- CM_plsda_training$confusionMatrices_mean
conf_mean_long <- melt(conf_mean, varnames = c("Reference", "Prediction"), value.name = "Mean")

# remove 0 values
conf_mean_long <- conf_mean_long[conf_mean_long$Mean > 0, ] 

# Define the color palette from white to orange
cols <- colorRampPalette(c('white', '#fe9929'))

# Plot using ggplot with Mean as the fill
cm_plot <- ggplot(data = conf_mean_long, aes(x = Prediction, y = Reference, fill = Mean)) +
  geom_tile(color = "white") +  # White border for tiles
  geom_text(aes(label = round(Mean, 1)), color = '#542788', size = 4) +  # Purple text for Mean values
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
  ggtitle("LDA Confusion Matrix for Herbarium Leaf Spectra") +
  labs(y = "True Identity", x = "Predicted Identity")

cm_plot

# Save the plot
ggsave("CM_LDA_ConfusionMatrix.pdf", plot = cm_plot, width = 7.5, height = 7.5)



#######
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
