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

root_path <- getwd()


#-------------------------------------------------------------------------------
#' @Spectra-Plots-Fig-S3
#-------------------------------------------------------------------------------

# plot of mean, quartiles (50%, 90%), and coeff. var. for each species

# Step 4: Create the plot
plot <- ggplot(statistics_data, aes(x = wavelength)) +
  # Add ribbons for distributions
  geom_ribbon(aes(ymin = low_90, ymax = high_90, fill = "90% Distribution"), alpha = 0.5) +
  geom_ribbon(aes(ymin = low_50, ymax = high_50, fill = "50% Distribution"), alpha = 0.9) +
  # Add lines for mean and CV
  geom_line(aes(y = mean, color = "Mean Reflectance"), size = 0.8) +
  geom_line(aes(y = cv, color = "Coefficient of Variation (CV)"), size = 0.8) +
  # Facet by species
  facet_wrap(~ species, scales = "free_y", ncol = 4) +
  # Labels and theme
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance / Coefficient of Variation",
    color = "Lines",  # Legend title for lines
    fill = "Shading"  # Legend title for ribbons
  ) +
  scale_fill_manual(values = c(
    "90% Distribution" = "gray80",
    "50% Distribution" = "gray50"
  )) +
  scale_color_manual(values = c(
    "Mean Reflectance" = "black",
    "Coefficient of Variation (CV)" = "red"
  )) +
  guides(
    color = guide_legend(order = 1),  # Place lines first
    fill = guide_legend(order = 2)   # Place shading second
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Customize axis text
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),  # Axis title size
    strip.text = element_text(size = 7),   # Facet label size
    legend.position = "bottom",           # Place legend at the bottom
    legend.box = "horizontal",            # Arrange legend items horizontally
    legend.title = element_blank(),       # Remove legend title
    legend.text = element_text(size = 8)  # Legend text size
  )

# Print and save the plot
print(plot)
ggsave("reflectance_statistics_by_species_reordered_legend.pdf", plot = plot, width = 12, height = 8)


####

# Step 1: Read the CSV file
specdf <- read.csv("data/dataHUH2024_sp25leaf561_ref5nm_450-2400.csv")

# Step 2: Melt the data into long format
melted_data <- specdf %>%
  select(name, species, 23:ncol(specdf)) %>%
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

# Step 4: Create the plot
specplot <- ggplot(statistics_data, aes(x = wavelength)) +
  # Add ribbons for distributions
  geom_ribbon(aes(ymin = low_90, ymax = high_90, fill = "90% Distribution")) +
  geom_ribbon(aes(ymin = low_50, ymax = high_50, fill = "50% Distribution")) +
  # Add lines for mean and CV
  geom_line(aes(y = mean, color = "Mean Reflectance"), size = 0.8) +
  geom_line(aes(y = cv, color = "Coefficient of Variation (CV)"), size = 0.8) +
  # Facet by species
  facet_wrap(~ species, scales = "free_y", ncol = 4) +
  # Labels and theme
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance / Coefficient of Variation",
    color = "Lines",  # Legend title for lines
    fill = "Shading"  # Legend title for ribbons
  ) +
  scale_fill_manual(values = c(
    "90% Distribution" = "gray80",
    "50% Distribution" = "gray50"
  )) +
  scale_color_manual(values = c(
    "Mean Reflectance" = "black",
    "Coefficient of Variation (CV)" = "red"
  )) +
  guides(
    color = guide_legend(order = 1),  # Place lines first
    fill = guide_legend(order = 2)   # Place shading second
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Customize axis text
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),  # Axis title size
    strip.text = element_text(size = 7),   # Facet label size
    legend.position = "bottom",           # Place legend at the bottom
    legend.box = "horizontal",            # Arrange legend items horizontally
    legend.title = element_blank(),       # Remove legend title
    legend.text = element_text(size = 8)  # Legend text size
  )

# Print and save the plot
print(specplot)
ggsave("Figures_Tables/Fig-S2_Spectra-CV.png", plot = specplot, width = 7, height = 8.5, dpi = 300)

#'------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------
