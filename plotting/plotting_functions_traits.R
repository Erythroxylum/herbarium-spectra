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

spectra_path <- "/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/herbarium-spectra"

#-------------------------------------------------------------------------------
#' @Trait-Biplots
#-------------------------------------------------------------------------------

###############
# Plot trait predictions

setwd("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/results/leaf_kothari5nm/")

# read obs-pred-avg huh
obs_pred_data <- fread("pls_huh5nm_testing_obs-pred-avg.csv")
obs_pred_data <- fread("pls_huh5nm_training_obs-pred-avg.csv")

# read obs-pred-avg Kothari
obs_pred_data <- fread("pls_kothari5nm_testing_obs-pred-avg.csv")
obs_pred_data <- fread("pls_kothari5nm_training_obs-pred-avg.csv")

# cross testing
obs_pred_data <- fread("pls_herbMod_presSpec_obs-pred-avg.csv")
obs_pred_data <- fread("pls_presMod_herbSpec_obs-pred-avg.csv")

plotname <- "pls_huh5nm_testing_obs-pred-avg"

# Create the plot
p <- ggplot(obs_pred_data, aes(x = predicted_mean, y = observed)) +
  geom_point(aes(color = species), alpha = 0.6) + # Points with some transparency
  geom_errorbarh(aes(xmin = predicted_mean - predicted_sd, 
                     xmax = predicted_mean + predicted_sd, 
                     color = species), # Match error bar color to species
                 height = 0.0) + # Horizontal error bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + # 1:1 line
  labs(
    title = plotname,
    x = expression("Predicted LMA (kg m"^-2*")"), # Superscript for -2
    y = expression("Measured LMA (kg m"^-2*")")
  ) +
  theme_minimal() + # Clean theme
  theme(legend.position = "right") +# Legend on the right
  guides(color = guide_legend(ncol = 1)) # Make the legend a single column

# Save the plot to a PDF file
ggsave(paste0(plotname,".pdf",sep=""), plot = p, width = 10, height = 8)


#-------------------------------------------------------------------------------
#' @Trait-Performance-Tables
#-------------------------------------------------------------------------------

# Load data
perf <- fread("pls_huh5nm_training_performance.csv")
perf <- fread("pls_huh5nm_testing_performance.csv")
perf <- fread("pls_kothari5nm_testing_performance.csv")
perf <- fread("pls_kothari5nm_training_performance.csv")
perf <- fread("pls_herbMod_presSpec_performance.csv")
perf <- fread("pls_presMod_herbSpec_performance.csv")

# Convert to data frame
perf_df <- as.data.frame(perf)

# Calculate mean and standard deviation for selected metrics
summary_stats <- perf_df %>%
  summarise(
    Nmodels = length(iteration),
    Nsamples = length(nsamp),
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2 = sd(R2, na.rm = TRUE),
    mean_BIAS = mean(BIAS, na.rm = TRUE),
    sd_BIAS = sd(BIAS, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_RMSE = sd(RMSE, na.rm = TRUE),
    mean_perRMSE = mean(perRMSE, na.rm = TRUE),
    sd_perRMSE = sd(perRMSE, na.rm = TRUE),
    mean_intercept = mean(intercept, na.rm = TRUE),
    sd_intercept = sd(intercept, na.rm = TRUE),
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE),
  )

# add split, models, and spectra info
summary_stats <- summary_stats %>%
  mutate(
    split = "All",
    models = "Pressed",
    spectra = "Herbarium",
  )

#Done
training_herbMod_herbSpec <- summary_stats
testing_herbMod_herbSpec <- summary_stats
testing_presMod_presMod <- summary_stats
training_presMod_presSpec <- summary_stats
all_herbMod_presSpec <- summary_stats
all_presMod_herbSpec <- summary_stats

#Pending


# Combine Stats for single table
combined_stats <- bind_rows(
  training_herbMod_herbSpec,
  testing_herbMod_herbSpec,
  training_presMod_presSpec,
  testing_presMod_presMod,
  all_presMod_herbSpec,
  all_herbMod_presSpec,
)

# Select and arrange columns in the specified order
final_table <- combined_stats %>%
  select(spectra, models, split, Nmodels, Nsamples,mean_R2, sd_R2, 
         mean_BIAS, sd_BIAS, mean_RMSE, sd_RMSE, mean_perRMSE, sd_perRMSE,
         mean_intercept, sd_intercept, mean_slope, sd_slope)

write.csv(final_table, "pls_performance_stats.csv")



#-------------------------------------------------------------------------------
#' @ PLS VIP Plots
#-------------------------------------------------------------------------------

library(data.table)

vip <- fread(paste0(root_path,"/pls_huh5nm_vip.csv"))
vip <- fread(paste0(root_path,"/pls_kothari5nm_vip.csv"))

vipname <- "pls_kothari5nm_vip"

# 'vip' is data table with wavelengths as columns and iterations as rows
# Remove the iteration column to calculate stats for wavelength columns
vip_wavelengths <- vip[, -1, with = FALSE]

# Calculate mean and coefficient of variation for each wavelength
mean_values <- colMeans(vip_wavelengths, na.rm = TRUE)
cv_values <- apply(vip_wavelengths, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))

# Create a data table for plotting with proper numeric conversion
vip_stats <- data.table(
  wavelength = as.numeric(gsub("[^0-9.]", "", names(mean_values))), # Remove non-numeric characters
  mean_VIP = mean_values,
  CV_VIP = cv_values
)

# Function to calculate desired percentiles
calculate_percentiles <- function(x) {
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
}

# Calculate the percentiles for each wavelength (skipping the iteration column)
percentile_data <- vip[, lapply(.SD, calculate_percentiles), .SDcols = -1]

# Convert this into a more usable format
# Create a list of quantile names (5%, 25%, 50%, 75%, 95%)
quantile_names <- c("low_90", "low_50", "mean", "high_50", "high_90")

# Transpose to long format
percentile_data_long <- as.data.table(t(percentile_data))
setnames(percentile_data_long, names(percentile_data_long), quantile_names)

# Add the corresponding wavelengths as a column
percentile_data_long[, wavelength := as.numeric(gsub("[^0-9.]", "", names(vip)[-1]))]  # Skip 'iteration' column

# Merge with the existing vip_stats (mean_VIP and CV_VIP) for plotting
vip_stats <- merge(vip_stats, percentile_data_long, by = "wavelength", all = TRUE)

# Check the structure of the final data
str(vip_stats)

library(ggplot2)

# Create the plot
p1 <- ggplot(vip_stats, aes(x = wavelength)) +
  geom_ribbon(aes(ymin = low_90, ymax = high_90), fill = "blue", alpha = 0.5) +  # 90% interval (shaded)
  geom_ribbon(aes(ymin = low_50, ymax = high_50), fill = "blue", alpha = 0.9) +  # 50% interval (shaded)
  geom_line(aes(y = mean_VIP, color = "Mean VIP"), linewidth = 0.8) +  # Mean VIP line
  geom_line(aes(y = CV_VIP, color = "CV VIP"), linewidth = 0.8) +  # CV VIP line
  labs(
    title = vipname,
    x = "Wavelength (nm)",
    y = "VIP Values",
    color = "Metric"  # Legend title
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.5),  # Major gridlines
    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25),  # Minor gridlines
    legend.position = "right"
  ) +
  scale_x_continuous(breaks = seq(0, max(vip_stats$wavelength, na.rm = TRUE), by = 200)) +  # Grid every 200 nm
  scale_y_continuous() +  # Default for y-axis grid
  scale_color_manual(values = c("Mean VIP" = "black", "CV VIP" = "red"))  # Manual color scale for lines

# Save the plot
ggsave(paste0(vipname,".pdf",sep=""), plot = p1, width = 8, height = 6)


