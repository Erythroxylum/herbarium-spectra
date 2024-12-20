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
library(patchwork) 

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

spectra_path <- getwd()

#-------------------------------------------------------------------------------
#' @Trait-Biplots-Fig1
#-------------------------------------------------------------------------------

##################
### HUH Data ref full (or refnorm full)
huh <- fread("results/HUH/ref/HUH_450-2400/pls_HUH_testing_obs-pred.csv")

# Group by 'accession_leaf' for HUH
averaged_data_huh <- huh %>%
  group_by(accession_leaf) %>%
  summarise(
    # Summarize all numerical columns with their mean
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    # Select the first value for character/logical columns
    across(where(is.character), first),
    across(where(is.logical), first)
  )

##################
### Kothari Data ref full (or refnorm full)
kothari <- fread("results/Kothari/ref/Kothari_450-2400/pls_Kothari_training_obs-pred.csv")

# Filter Kothari to include only species present in the HUH dataset
huh_species <- unique(huh$species)
kothari_filtered <- kothari %>% filter(species %in% huh_species)

# Group by 'accession' for Kothari
averaged_data_kothari <- kothari_filtered %>%
  group_by(accession) %>%
  summarise(
    # Summarize all numerical columns with their mean
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    # Select the first value for character/logical columns
    across(where(is.character), first),
    across(where(is.logical), first)
  )


##################
### Model Herb Spec Press cwt 450-2400 (or ref full)
mhsp <- fread("results/model_transfer/LMA_coef-HUH_spec-Kothari_cwt_450-2400/coef-HUH_spec-Kothari_cwt_450-2400_obs-pred.csv")

# Filter mpsh to include only species present in the HUH dataset
huh_species <- unique(huh$species)
mhsp_filtered <- mhsp %>% filter(species %in% huh_species)

# Group by 'accession' for mpsh
averaged_data_mhsp <- mhsp_filtered %>%
  group_by(accession) %>%
  summarise(
    # Summarize all numerical columns with their mean
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    # Select the first value for character/logical columns
    across(where(is.character), first),
    across(where(is.logical), first)
  )


##################
### Model press Spec Herb cwt 450-2400 (or ref 450-1300)
mpsh <- fread("results/model_transfer/LMA_coef-Kothari_spec-HUH_cwt_450-2400/coef-Kothari_spec-HUH_cwt_450-2400_obs-pred.csv")

# Filter mpsh to include only species present in the HUH dataset
huh_species <- unique(huh$species)
mpsh_filtered <- mpsh %>% filter(species %in% huh_species)

# Group by 'accession' for mpsh
averaged_data_mpsh <- mpsh_filtered %>%
  group_by(accession) %>%
  summarise(
    # Summarize all numerical columns with their mean
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    # Select the first value for character/logical columns
    across(where(is.character), first),
    across(where(is.logical), first)
  )


##################
### Define Consistent Colors
# Combine all unique species across datasets
all_species <- sort(unique(c(huh$species, kothari_filtered$species)))

# Generate consistent colors for all species
species_colors <- setNames(viridisLite::turbo(length(all_species)), all_species)


##################
### Create a light gray theme for all plots
custom_theme <- theme(
  panel.background = element_rect(fill = "gray93", color = NA),
  legend.position = "none", # No legend
  axis.text = element_text(size = 8),  # Uniform axis text size
  axis.title = element_text(size = 8) # Uniform axis title size
)

##################
### Kothari Plot
kothari_plot <- ggplot(averaged_data_kothari, aes(x = predicted_mean, y = observed)) +
  geom_point(aes(color = species), alpha = 0.6) +
  geom_errorbarh(aes(xmin = predicted_mean - predicted_sd, xmax = predicted_mean + predicted_sd, color = species), height = 0.0) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = NULL, y = expression("Measured LMA (kg m"^-2*")")) +
  scale_color_manual(values = species_colors) +
  custom_theme +
  coord_fixed(ratio = 1, xlim = c(0, 0.225), ylim = c(0, 0.225)) +
  annotate("text", x = 0.01, y = 0.21, label = "Test: Validation\nModel: pressed\nSpectra: pressed\nTransform: reflectance", size = 2.5, hjust = 0)


##################
### HUH Plot with legend
huh_plot <- ggplot(averaged_data_huh, aes(x = predicted_mean, y = observed)) +
  geom_point(aes(color = species), alpha = 0.6) +
  geom_errorbarh(aes(xmin = predicted_mean - predicted_sd, xmax = predicted_mean + predicted_sd, color = species), height = 0.0) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = NULL, y = NULL) +
  scale_color_manual(values = species_colors) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.15, "cm"),
    legend.text = element_text(size = 8, face = "italic"),
    legend.title = element_text(size = 8),
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm"),
    panel.background = element_rect(fill = "gray93", color = NA),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  ) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  coord_fixed(ratio = 1, xlim = c(0, 0.225), ylim = c(0, 0.225)) +
  annotate("text", x = 0.01, y = 0.21, label = "Test: Validation\nModel: herbarium\nSpectra: herbarium\nTransform: reflectance", size = 2.5, hjust = 0)

##################
### MPSH Plot
mpsh_plot <- ggplot(averaged_data_mpsh, aes(x = predicted_mean, y = observed)) +
  geom_point(aes(color = species), alpha = 0.6) +
  geom_errorbarh(aes(xmin = predicted_mean - predicted_sd, xmax = predicted_mean + predicted_sd, color = species), height = 0.0) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = expression("Predicted LMA (kg m"^-2*")"), y = NULL) +
  scale_color_manual(values = species_colors) +
  custom_theme +
  coord_fixed(ratio = 1, xlim = c(0, 0.225), ylim = c(0, 0.225)) +
  annotate("text", x = 0.01, y = 0.21, label = "Test: Transfer\nModel: pressed\nSpectra: herbarium\nTransform: CWT", size = 2.5, hjust = 0)


##################
### MHSP Plot
mhsp_plot <- ggplot(averaged_data_mhsp, aes(x = predicted_mean, y = observed)) +
  geom_point(aes(color = species), alpha = 0.6) +
  geom_errorbarh(aes(xmin = predicted_mean - predicted_sd, xmax = predicted_mean + predicted_sd, color = species), height = 0.0) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = expression("Predicted LMA (kg m"^-2*")"), y = expression("Measured LMA (kg m"^-2*")")) +
  scale_color_manual(values = species_colors) +
  custom_theme +
  coord_fixed(ratio = 1, xlim = c(0, 0.225), ylim = c(0, 0.225)) +
  annotate("text", x = 0.01, y = 0.21, label = "Test: Transfer\nModel: herbarium\nSpectra: pressed\nTransform: CWT", size = 2.5, hjust = 0)


##################
### Combine Plots 2 x 2
combined_plot <- wrap_plots(
  kothari_plot, huh_plot, mhsp_plot, mpsh_plot, 
  ncol = 2, nrow = 2, guides = "collect"
) +
  plot_annotation(tag_levels = 'A')

# Save the combined plot to a PDF
ggsave("results/Fig2_LMA-biplots.pdf", plot = combined_plot, width = 8.5, height = 6.5)







#-------------------------------------------------------------------------------
#' @PLS-VIP-Fig
#-------------------------------------------------------------------------------
library(data.table)
library(ggplot2)
library(patchwork)

root_path <- getwd()

# Function to read and process VIP data
process_vip <- function(file_path) {
  # Read the VIP data
  vip <- fread(file_path)
  
  # Remove the iteration column to calculate stats for wavelength columns
  vip_wavelengths <- vip[, -1, with = FALSE]
  
  # Calculate mean and coefficient of variation for each wavelength
  mean_values <- colMeans(vip_wavelengths, na.rm = TRUE)
  cv_values <- apply(vip_wavelengths, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  
  # Create a data table for plotting with proper numeric conversion
  vip_stats <- data.table(
    wavelength = as.numeric(gsub("[^0-9.]", "", names(mean_values))), # Extract numeric wavelengths
    mean_VIP = mean_values,
    CV_VIP = cv_values
  )
  
  # Function to calculate desired percentiles
  calculate_percentiles <- function(x) {
    quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
  }
  
  # Calculate the percentiles for each wavelength (skipping the iteration column)
  percentile_data <- vip[, lapply(.SD, calculate_percentiles), .SDcols = -1]
  
  # Create a list of quantile names (5%, 25%, 50%, 75%, 95%)
  quantile_names <- c("low_90", "low_50", "mean", "high_50", "high_90")
  
  # Transpose to long format
  percentile_data_long <- as.data.table(t(percentile_data))
  setnames(percentile_data_long, names(percentile_data_long), quantile_names)
  
  # Add the corresponding wavelengths as a column
  percentile_data_long[, wavelength := as.numeric(gsub("[^0-9.]", "", names(vip)[-1]))]  # Skip 'iteration' column
  
  # Merge with the existing vip_stats (mean_VIP and CV_VIP) for plotting
  vip_stats <- merge(vip_stats, percentile_data_long, by = "wavelength", all = TRUE)
  
  return(vip_stats)
}

# Function to create a plot from processed VIP data
plot_vip <- function(vip_stats, title) {
  ggplot(vip_stats, aes(x = wavelength)) +
    geom_ribbon(aes(ymin = low_90, ymax = high_90), fill = "blue", alpha = 0.5) +  # 90% interval (shaded)
    geom_ribbon(aes(ymin = low_50, ymax = high_50), fill = "blue", alpha = 0.9) +  # 50% interval (shaded)
    geom_line(aes(y = mean_VIP, color = "Mean VIP"), linewidth = 0.8) +  # Mean VIP line
    geom_line(aes(y = CV_VIP, color = "CV VIP"), linewidth = 0.8) +  # CV VIP line
    labs(
      title = title,
      x = "Wavelength (nm)",
      y = "VIP Values",
      color = "Metric"  # Legend title
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray", linewidth = 0.5),  # Major gridlines
      panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25),  # Minor gridlines
      legend.position = "none", # or "right"
      plot.title = element_text(size = 12, hjust = 0.1)
    ) +
    scale_x_continuous(breaks = seq(0, max(vip_stats$wavelength, na.rm = TRUE), by = 200)) +  # Grid every 200 nm
    scale_y_continuous() +
    scale_color_manual(values = c("Mean VIP" = "black", "CV VIP" = "red"))  # Manual color scale for lines
}

# List of datasets, output names, and titles
datasets <- list(
  HUHref = paste0(root_path, "/results/HUH/ref/HUH_450-2400/pls_HUH_vip.csv"),
  Kotref = paste0(root_path, "/results/Kothari/ref/Kothari_450-2400/pls_Kothari_vip.csv"),
  HUHcwt = paste0(root_path, "/results/HUH/cwt/HUH_450-2400/pls_HUH_vip.csv"),
  Kotcwt = paste0(root_path, "/results/Kothari/cwt/Kothari_450-2400/pls_Kothari_vip.csv")
)

titles <- list(
  HUHref = "Model and spectra: herbarium, Transform: reflectance",
  Kotref = "Model and spectra: pressed, Transform: reflectance",
  HUHcwt = "Model and spectra: herbarium, Transform: CWT",
  Kotcwt = "Model and spectra: pressed, Transform: CWT"
)

# Process and plot each dataset with corresponding titles
plots <- lapply(names(datasets), function(name) {
  vip_stats <- process_vip(datasets[[name]])
  plot_vip(vip_stats, title = titles[[name]])
})

# Combine plots into a single figure
vip_plot <- wrap_plots(plots, ncol = 1, guides = "collect") +
  plot_annotation(tag_levels = 'A')

# Save the combined plot
ggsave("results/Fig3_VIP.pdf", plot = vip_plot, width = 8.5, height = 9)


#-------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------

