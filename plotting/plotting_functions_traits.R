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

root_path <- getwd()

#-------------------------------------------------------------------------------
#' @PLS-VIP-sixPanel-Fig3
#-------------------------------------------------------------------------------

## Process Datasets

process_reflectance <- function(file_path, bands, band_type) {
  # Read the full dataset
  full_data <- fread(file_path)
  
  # Ensure bands match column names
  numeric_bands <- as.character(bands)  # Convert bands to character
  
  # Select only reflectance columns
  reflectance_data <- full_data[, ..numeric_bands, with = FALSE]
  
  # Calculate mean and percentiles for each wavelength
  mean_values <- colMeans(reflectance_data, na.rm = TRUE)
  percentile_values <- apply(reflectance_data, 2, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))
  
  # Create a data table for plotting
  data.table(
    wavelength = as.numeric(names(mean_values)),  # Wavelengths are column names
    mean = mean_values,
    low_90 = percentile_values[1, ],
    high_90 = percentile_values[2, ],
    type = band_type  # Add band type (Reflectance or CWT Reflectance)
  )
}


process_vip <- function(file_path, band_type) {
  # Read the dataset
  vip_coeff_data <- fread(file_path)
  
  # Remove "iteration" and "(Intercept)" columns
  vip_coeff_data <- vip_coeff_data[, !c("iteration"), with = FALSE]
  
  # Clean column names (remove backticks and ensure numeric)
  setnames(vip_coeff_data, gsub("`", "", colnames(vip_coeff_data)))
  
  # Calculate mean and percentiles for each wavelength
  mean_values <- colMeans(vip_coeff_data, na.rm = TRUE)
  percentile_values <- apply(vip_coeff_data, 2, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))
  
  # Create a data table for plotting
  data.table(
    wavelength = as.numeric(names(mean_values)),  # Convert cleaned column names to numeric
    mean = mean_values,
    low_90 = percentile_values[1, ],
    high_90 = percentile_values[2, ],
    type = band_type  # Add band type (VIP or Coefficients)
  )
}

process_coeff <- function(file_path, band_type) {
  # Read the dataset
  vip_coeff_data <- fread(file_path)
  
  # Remove "iteration" and "(Intercept)" columns
  vip_coeff_data <- vip_coeff_data[, !c("iteration", "(Intercept)"), with = FALSE]
  
  # Clean column names (remove backticks and ensure numeric)
  setnames(vip_coeff_data, gsub("`", "", colnames(vip_coeff_data)))
  
  # Calculate mean and percentiles for each wavelength
  mean_values <- colMeans(vip_coeff_data, na.rm = TRUE)
  percentile_values <- apply(vip_coeff_data, 2, function(x) quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))
  
  # Create a data table for plotting
  data.table(
    wavelength = as.numeric(names(mean_values)),  # Convert cleaned column names to numeric
    mean = mean_values,
    low_90 = percentile_values[1, ],
    high_90 = percentile_values[2, ],
    type = band_type  # Add band type (VIP or Coefficients)
  )
}



## Read in datasets

# Define bands for reflectance datasets
#bands <- paste0("X", seq(450, 2400, by = 5)) 
bands <- seq(450, 2400, by = 5) # Match column naming convention

reflectance_data <- rbindlist(list(
  process_reflectance(paste0(root_path, "/data/dataHUH2024_sp25leaf561_ref5nm_450-2400.csv"), bands, "Herbarium"),
  process_reflectance(paste0(root_path, "/data/dataKothari_pressed_unavg_ref5nm_450-2400.csv"), bands, "Pressed")
))

# CWT Reflectance Data
cwt_reflectance_data <- rbindlist(list(
  process_reflectance(paste0(root_path, "/data/dataHUH2024_sp25leaf561_cwt5nm_450-2400.csv"), bands, "Herbarium"),
  process_reflectance(paste0(root_path, "/data/dataKothari_pressed_unavg_cwt5nm_450-2400.csv"), bands, "Pressed")
))

# VIP Data
vip_data <- rbindlist(list(
  process_vip(paste0(root_path, "/results/HUH/ref/HUH_450-2400/pls_HUH_vip.csv"), "Herbarium"),
  process_vip(paste0(root_path, "/results/Kothari/ref/Kothari_450-2400/pls_Kothari_vip.csv"), "Pressed")
))

# Coefficients Data
coeff_data <- rbindlist(list(
  process_coeff(paste0(root_path, "/results/HUH/ref/HUH_450-2400/pls_HUH_coefficients.csv"), "Herbarium"),
  process_coeff(paste0(root_path, "/results/Kothari/ref/Kothari_450-2400/pls_Kothari_coefficients.csv"), "Pressed")
))

# Process CWT VIP data
cwt_vip_data <- rbindlist(list(
  process_vip(
    paste0(root_path, "/results/HUH/cwt/HUH_450-2400/pls_HUH_vip.csv"),
    band_type = "Herbarium"
  ),
  process_vip(
    paste0(root_path, "/results/Kothari/cwt/Kothari_450-2400/pls_Kothari_vip.csv"),
    band_type = "Pressed"
  )
))

# CWT Coefficients Data
cwt_coeff_data <- rbindlist(list(
  process_coeff(paste0(root_path, "/results/HUH/cwt/HUH_450-2400/pls_HUH_coefficients.csv"), "Herbarium"),
  process_coeff(paste0(root_path, "/results/Kothari/cwt/Kothari_450-2400/pls_Kothari_coefficients.csv"), "Pressed")
))

## Plotting Function

create_panel_plot <- function(data, color_map, x_title = TRUE, y_title = TRUE) {
  ggplot(data, aes(x = wavelength, group = type)) +
    geom_ribbon(aes(ymin = low_90, ymax = high_90, fill = "90% quantile"), alpha = 0.5, color = NA) +  # Gray band
    geom_line(aes(y = mean, color = type), size = 0.5) +  # Mean line
    scale_color_manual(
      values = color_map,
      labels = c("Herbarium mean value", "Pressed mean value")  # Custom legend labels for lines
    ) +
    scale_fill_manual(
      values = c("90% quantile" = "gray70"),
      labels = c("90% quantile")  # Custom legend label for gray ribbon
    ) +
    guides(
      color = guide_legend(order = 1),  # Place lines first in the legend
      fill = guide_legend(order = 2)    # Place ribbon second in the legend
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",  # Disable legend for individual plots
      axis.title.x = if (x_title) element_text(size = 10) else element_blank(),  # Conditionally set x-axis title
      axis.text.x = if (x_title) element_text(size = 8) else element_blank(),  # Conditionally set x-axis tick labels
      axis.title.y = if (y_title) element_text(size = 10) else element_blank(),  # Conditionally set y-axis title
      axis.text.y = element_text(size = 8),  # Always display y-axis tick labels
      strip.text = element_text(size = 10),
      plot.title = element_blank()  # Always remove plot titles
    ) +
    labs(x = if (x_title) "Wavelength (nm)" else NULL, y = if (y_title) "Value" else NULL)  # Dynamically set axis labels
}


## create panels

# Panel A and B (Reflectance)
panel_a <- create_panel_plot(reflectance_data, c("Herbarium" = "black", "Pressed" = "red"), x_title = FALSE, y_title = TRUE)
panel_b <- create_panel_plot(cwt_reflectance_data, c("Herbarium" = "black", "Pressed" = "red"), x_title = FALSE, y_title = FALSE)

# Panel C and D (VIP)
panel_c <- create_panel_plot(vip_data, c("Herbarium" = "black", "Pressed" = "red"), x_title = FALSE, y_title = TRUE)
panel_d <- create_panel_plot(cwt_vip_data, c("Herbarium" = "black", "Pressed" = "red"), x_title = FALSE, y_title = FALSE)

# Panel E and F (Coefficients)
panel_e <- create_panel_plot(coeff_data, c("Herbarium" = "black", "Pressed" = "red"), x_title = TRUE, y_title = TRUE)
panel_f <- create_panel_plot(cwt_coeff_data, c("Herbarium" = "black", "Pressed" = "red"), x_title = TRUE, y_title = FALSE)


# Combine plots into a single figure
final_figure <- wrap_plots(
  panel_a, panel_b,
  panel_c, panel_d,
  panel_e, panel_f,
  ncol = 2,  # Arrange into 2 columns
  guides = "collect"  # Collect legends from all panels
) +
  plot_annotation(tag_levels = "A") &  # Add annotations (A, B, C, ...)
  theme(
    legend.position = "bottom",  # Set legend position at the bottom
    legend.box = "horizontal",  # Arrange legend items horizontally
    legend.title = element_blank()  # Remove legend title
  )

# Save the final figure
ggsave("Figures_Tables/Fig3-6panel.pdf", plot = final_figure, width = 8.5, height = 6)

#-------------------------------------------------------------------------------
#' @Trait-Biplots-Fig4
#-------------------------------------------------------------------------------

##################
### HUH Data ref full (or refnorm full)
huh <- fread("results/HUH/ref/HUH_450-2400/pls_HUH_testing_obs-pred.csv")

# Group by 'accession_leaf' for HUH
averaged_data_huh <- huh %>%
  group_by(accession) %>%
  summarise(
    # Summarize all numerical columns with their mean
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    # Select the first value for character/logical columns
    across(where(is.character), first),
    across(where(is.logical), first)
  )

##################
### Kothari Data ref full (or refnorm full)
kothari <- fread("results/Kothari/ref/Kothari_450-2400/pls_Kothari_testing_obs-pred.csv")

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

mhsp_line <- fread("results/model_transfer/LMA_coef-HUH_spec-Kothari_cwt_450-2400/coef-HUH_spec-Kothari_cwt_450-2400_performance.csv")
slope <- mean(mhsp_line$slope)
intercept <- mean(mhsp_line$intercept)

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


# Kothari Plot
kothari_plot <- ggplot(averaged_data_kothari, aes(x = predicted_mean, y = observed)) +
  geom_point(aes(color = species), alpha = 0.6) +
  geom_errorbarh(aes(xmin = predicted_mean - predicted_sd, xmax = predicted_mean + predicted_sd, color = species), height = 0.0) +
  #geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") +  # Regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray", size=1.5) +  # 1:1 line
  geom_abline(slope = 1.017387475, intercept = -0.00323509, linetype = "solid", color = "red", size = 0.5) +  # Custom line
  labs(x = NULL, y = expression("Measured LMA (kg m"^-2*")")) +
  scale_color_manual(values = species_colors) +
  custom_theme +
  coord_fixed(ratio = 1, xlim = c(0, 0.225), ylim = c(0, 0.225)) +
  annotate("text", x = 0.01, y = 0.21, label = "Test: validation\nModel: pressed\nSpectra: pressed\nTransform: reflectance", size = 2.5, hjust = 0)

# HUH Plot with legend
huh_plot <- ggplot(averaged_data_huh, aes(x = predicted_mean, y = observed)) +
  geom_point(aes(color = species), alpha = 0.6) +
  geom_errorbarh(aes(xmin = predicted_mean - predicted_sd, xmax = predicted_mean + predicted_sd, color = species), height = 0.0) +
  #geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") +  # Regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray", size = 1.5) +  # 1:1 line
  geom_abline(slope = 0.99, intercept = 7.88e-04, linetype = "solid", color = "red", size = 0.5) +  # Custom line
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
  annotate("text", x = 0.01, y = 0.21, label = "Test: validation\nModel: herbarium\nSpectra: herbarium\nTransform: reflectance", size = 2.5, hjust = 0)

# MPSH Plot
mpsh_plot <- ggplot(averaged_data_mpsh, aes(x = predicted_mean, y = observed)) +
  geom_point(aes(color = species), alpha = 0.6) +
  geom_errorbarh(aes(xmin = predicted_mean - predicted_sd, xmax = predicted_mean + predicted_sd, color = species), height = 0.0) +
  #geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") +  # Regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray", size = 1.5) +  # 1:1 line
  geom_abline(slope = 1.252121176, intercept = -0.021495026, linetype = "solid", color = "red", size = 0.5) +  # Custom line
  labs(x = expression("Predicted LMA (kg m"^-2*")"), y = NULL) +
  scale_color_manual(values = species_colors) +
  custom_theme +
  coord_fixed(ratio = 1, xlim = c(0, 0.225), ylim = c(0, 0.225)) +
  annotate("text", x = 0.01, y = 0.21, label = "Test: transfer\nModel: pressed\nSpectra: herbarium\nTransform: CWT", size = 2.5, hjust = 0)

# MHSP Plot
# intercept:0.004336495, slope:0.9142096
mhsp_plot <- ggplot(averaged_data_mhsp, aes(x = predicted_mean, y = observed)) +
  geom_point(aes(color = species), alpha = 0.6) +
  geom_errorbarh(aes(xmin = predicted_mean - predicted_sd, xmax = predicted_mean + predicted_sd, color = species), height = 0.0) +
  #geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") +  # Regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray", size = 1.5) +  # 1:1 line
  geom_abline(slope = 0.9142096, intercept = 0.004336495, linetype = "solid", color = "red", size = 0.5) +  # Custom line
  labs(x = expression("Predicted LMA (kg m"^-2*")"), y = expression("Measured LMA (kg m"^-2*")")) +
  scale_color_manual(values = species_colors) +
  custom_theme +
  coord_fixed(ratio = 1, xlim = c(0, 0.225), ylim = c(0, 0.225)) +
  annotate("text", x = 0.01, y = 0.21, label = "Test: transfer\nModel: herbarium\nSpectra: pressed\nTransform: CWT", size = 2.5, hjust = 0)


##################
### Combine Plots 2 x 2
combined_plot <- wrap_plots(
  kothari_plot, huh_plot, mhsp_plot, mpsh_plot, 
  ncol = 2, nrow = 2, guides = "collect"
) +
  plot_annotation(tag_levels = 'A')

# Save the combined plot to a PDF
ggsave("Figures_Tables/Fig4_LMA-biplots.png", plot = combined_plot, width = 8.5, height = 6.5)


#-------------------------------------------------------------------------------
#' @Violin-plots-Fig5
#-------------------------------------------------------------------------------

# do cwt then ref then cwt
# do each HUH trait for each ggplot and then make final plot together.

# load Pressed data, any transformation.
Kothari_data <- read.csv("data/dataKothari_pressed_unavg_cwt5nm_450-2400.csv")

# set trait:  nothing or leafKg_m2 for LMA, C, Ca, car, cel or cellulose, chlA, N, sol or solubles
trait_name <- ""

# Step 1: Load HUH trait data - edit dataset "cwt 1350-2400" or "refnorm" or full range
HUH_data <- read.csv(paste0("results/model_transfer/Coef-Kothari_Spec-HUH_cwt", trait_name, "_450-1300/Coef-Kothari_Spec-HUH_cwt", trait_name, "_450-1300_predicted_trait.csv"))

# Step 2: Parse HUH Data
HUH_parse <- HUH_data %>%
  select(species, predicted_trait_mean) %>%
  mutate(dataset = factor("Herbarium", levels = c("Herbarium", "Pressed")), 
         trait_value = predicted_trait_mean) %>%
  select(species, trait_value, dataset)

# Step 3: Parse Kothari Trait Data, needs "leafKg_m2", "cellulose", or "solubles" instead
Kothari_filtered <- Kothari_data %>%
  filter(species %in% HUH_parse$species) %>%
  select(species, trait_value = all_of(trait_name)) %>%
  mutate(dataset = factor("Pressed", levels = c("Herbarium", "Pressed"))) %>% # Explicit factor levels
  select(species, trait_value, dataset)

# Combine the two datasets
combined_data_cwt_LMA <- bind_rows(HUH_parse, Kothari_filtered) %>%
  arrange(species, dataset)
combined_data_cwt_C <- bind_rows(HUH_parse, Kothari_filtered) %>%
  arrange(species, dataset)
combined_data_cwt_Ca <- bind_rows(HUH_parse, Kothari_filtered) %>%
  arrange(species, dataset)
combined_data_cwt_car <- bind_rows(HUH_parse, Kothari_filtered) %>%
  arrange(species, dataset)
combined_data_cwt_cel <- bind_rows(HUH_parse, Kothari_filtered) %>%
  arrange(species, dataset)
combined_data_cwt_chlA <- bind_rows(HUH_parse, Kothari_filtered) %>%
  arrange(species, dataset)
combined_data_cwt_N <- bind_rows(HUH_parse, Kothari_filtered) %>%
  arrange(species, dataset)
combined_data_cwt_sol <- bind_rows(HUH_parse, Kothari_filtered) %>%
  arrange(species, dataset)


## Step 4: Create a split violin plot for trait value distributions
# set combined data and scale_y_continuous name

# plots with no x axis labels
trait_LMA <- ggplot(combined_data_cwt_LMA, aes(x = species, y = trait_value, fill = dataset)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE, scale = "width") + # Adjust scale for wider violins
  stat_summary(fun.data = "mean_se", geom = "pointrange", 
               position = position_dodge(0.8), show.legend = FALSE, size = 0.1) + # Adjust size for visibility
  scale_x_discrete(name = NULL) +       # Remove x-axis title
  scale_y_continuous(name = expression(LMA ~ (Kg/m^2))) + # Superscript for y-axis label
  scale_fill_manual(values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      # Remove x-axis labels
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    legend.position = "none",           # Remove the legend
  )

trait_C <- ggplot(combined_data_cwt_C, aes(x = species, y = trait_value, fill = dataset)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE, scale = "width") + # Adjust scale for wider violins
  stat_summary(fun.data = "mean_se", geom = "pointrange", 
               position = position_dodge(0.8), show.legend = FALSE, size = 0.1) + # Adjust size for visibility
  scale_x_discrete(name = NULL) +       # Remove x-axis title
  scale_y_continuous(name = expression(C ~ (mg ~ g^{-1}))) + # Superscript for y-axis label
  scale_fill_manual(values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      # Remove x-axis labels
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    legend.position = "none",           # Remove the legend
  )

trait_Ca <- ggplot(combined_data_cwt_Ca, aes(x = species, y = trait_value, fill = dataset)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE, scale = "width") + # Adjust scale for wider violins
  stat_summary(fun.data = "mean_se", geom = "pointrange", 
               position = position_dodge(0.8), show.legend = FALSE, size = 0.1) + # Adjust size for visibility
  scale_x_discrete(name = NULL) +       # Remove x-axis title
  scale_y_continuous(name = expression(Ca ~ (mg ~ g^{-1}))) + # Superscript for y-axis label
  scale_fill_manual(values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      # Remove x-axis labels
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    legend.position = "none",           # Remove the legend
  )

trait_car <- ggplot(combined_data_cwt_car, aes(x = species, y = trait_value, fill = dataset)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE, scale = "width") + # Adjust scale for wider violins
  stat_summary(fun.data = "mean_se", geom = "pointrange", 
               position = position_dodge(0.8), show.legend = FALSE, size = 0.1) + # Adjust size for visibility
  scale_x_discrete(name = NULL) +       # Remove x-axis title
  scale_y_continuous(name = expression(Carotenoids ~ (mg ~ g^{-1}))) + # Superscript for y-axis label
  scale_fill_manual(values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      # Remove x-axis labels
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    legend.position = "none",           # Remove the legend
  )

trait_cel <- ggplot(combined_data_cwt_cel, aes(x = species, y = trait_value, fill = dataset)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE, scale = "width") + # Adjust scale for wider violins
  stat_summary(fun.data = "mean_se", geom = "pointrange", 
               position = position_dodge(0.8), show.legend = FALSE, size = 0.1) + # Adjust size for visibility
  scale_x_discrete(name = NULL) +       # Remove x-axis title
  scale_y_continuous(name = expression(Cellulose ~ (mg ~ g^{-1}))) + # Superscript for y-axis label
  scale_fill_manual(values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      # Remove x-axis labels
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7), 
    legend.position = "none",           # Remove the legend
  )

trait_chlA <- ggplot(combined_data_cwt_chlA, aes(x = species, y = trait_value, fill = dataset)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE, scale = "width") + # Adjust scale for wider violins
  stat_summary(fun.data = "mean_se", geom = "pointrange", 
               position = position_dodge(0.8), show.legend = FALSE, size = 0.1) + # Adjust size for visibility
  scale_x_discrete(name = NULL) +       # Remove x-axis title
  scale_y_continuous(name = expression(Chl ~ italic(a) ~ (mg ~ g^{-1}))) + # Superscript for y-axis label
  scale_fill_manual(values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      # Remove x-axis labels
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7), 
    legend.position = "none",           # Remove the legend
  )

trait_N <- ggplot(combined_data_cwt_N, aes(x = species, y = trait_value, fill = dataset)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE, scale = "width") + # Adjust scale for wider violins
  stat_summary(fun.data = "mean_se", geom = "pointrange", 
               position = position_dodge(0.8), show.legend = FALSE, size = 0.1) + # Adjust size for visibility
  scale_x_discrete(name = NULL) +       # Remove x-axis title
  scale_y_continuous(name = expression(N ~ (mg ~ g^{-1}))) + # Superscript for y-axis label
  scale_fill_manual(values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      # Remove x-axis labels
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 7), 
    legend.position = "none",           # Remove the legend
  )

## Bottom plot, with x axis legend
trait_sol <- ggplot(combined_data_cwt_sol, aes(x = species, y = trait_value, fill = dataset)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE, scale = "width") + # Adjust scale for wider violins
  stat_summary(fun.data = "mean_se", geom = "pointrange", 
               position = position_dodge(0.8), show.legend = FALSE, size = 0.1) + # Adjust size for visibility
  scale_x_discrete(labels = function(x) lapply(x, function(label) bquote(italic(.(label))))) + # Italicize species names
  scale_y_continuous(name = expression(Solubles ~ (mg ~ g^{-1}))) + # Correct y-axis for trait_value (continuous)
  scale_fill_manual(name = "Dataset: CWT, 450-1,300 nm ", values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8), # Format x-axis labels
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank(), # Remove x-axis title
    axis.title.y = element_text(size = 7), 
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )


##patchwork
combined_plot <- trait_LMA / trait_C / trait_Ca / trait_car / trait_cel / trait_chlA / trait_N / trait_sol

ggsave("results/Fig-4-trait_violins_cwt_1350-2400v2.pdf", plot = combined_plot, width = 8.5, height = 11)
ggsave("results/Fig-S-trait_violins_cwt_450-2400v2.pdf", plot = combined_plot, width = 8.5, height = 11)
ggsave("results/Fig-S-trait_violins_refnorm_1350-2400v2.pdf", plot = combined_plot, width = 8.5, height = 11)



#'------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------
