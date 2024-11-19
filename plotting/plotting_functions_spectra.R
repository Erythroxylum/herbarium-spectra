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

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# Select the root_folder to read data and export results
# Antonio
root_path <- "C:/Users/jog4076/Downloads"
# DW
root_path <- "/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/herbarium-spectra"


#-------------------------------------------------------------------------------
#' @Spectra-Plots
#-------------------------------------------------------------------------------

# plot of mean, quartiles (50%, 90%), and coeff. var. for each species

# Step 1: Read the CSV file
specdf <- read.csv("../fullDataHUH2024_sp25leaf636_noResample_400-2300_noVNorm.csv")

# Step 2: Melt the data into long format
melted_data <- specdf %>%
  select(name, species, 21:ncol(specdf)) %>%
  pivot_longer(cols = starts_with("X"), names_to = "wavelength", values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(sub("X", "", wavelength))) # Remove the 'X' and convert to numeric

# Step 3: Calculate statistics for each wavelength
statistics_data <- melted_data %>%
  group_by(species, wavelength) %>%
  summarise(
    mean = mean(reflectance, na.rm = TRUE),
    low_90 = quantile(reflectance, 0.05, na.rm = TRUE),
    low_50 = quantile(reflectance, 0.25, na.rm = TRUE),
    high_50 = quantile(reflectance, 0.75, na.rm = TRUE),
    high_90 = quantile(reflectance, 0.95, na.rm = TRUE),
    cv = (sd(reflectance, na.rm = TRUE) / mean(reflectance, na.rm = TRUE)),  # CV divided by 100
    .groups = "drop"
  )

# Step 4: Plot using ggplot2
ggplot(statistics_data, aes(x = wavelength)) +
  geom_ribbon(aes(ymin = low_50, ymax = high_50), fill = "grey", alpha = 0.9) + # Quartiles as a band (50th percentiles)
  geom_ribbon(aes(ymin = low_90, ymax = high_90), fill = "lightgrey", alpha = 0.5) + # Quartiles as a band (90th percentiles)
  geom_line(aes(y = cv), color = "red", size = 0.8) + # CV line
  geom_line(aes(y = mean), color = "black", size = 0.8) + # Mean line
  facet_wrap(~ species, scales = "free_y") + # Separate plots by species
  labs(
    #title = "Reflectance Statistics by Wavelength",
    x = "Wavelength (nm)",
    y = "Reflectance / Coefficient of Variation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 8), # Reduce size of axis text
    axis.title = element_text(size = 10), # Reduce size of axis titles
    strip.text = element_text(size = 7), # Reduce size of facet labels
    plot.title = element_text(size = 12, hjust = 0.5), # Reduce size of plot title
    plot.margin = margin(1, 1, 1, 1) # Ensure plot doesn't cut off content
  )

# Save the plot
ggsave("Fig_Spectra-CV.pdf", width = 8, height = 7)

