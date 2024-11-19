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
#' @ PLSDA VIP Plots
#-------------------------------------------------------------------------------

vip_plsda <- fread(paste0(root_path,"/plsda_varImp.csv"))

# Convert the data table to a data frame (as you've done)
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

# Step 4: Calculate mean, 50% and 90% quantiles, and CV for each Species at each Wavelength
vip_stats <- df_long[, .(
  Mean_VIP = mean(VIP_Value, na.rm = TRUE),
  low_50 = quantile(VIP_Value, 0.25, na.rm = TRUE),
  high_50 = quantile(VIP_Value, 0.75, na.rm = TRUE),
  low_90 = quantile(VIP_Value, 0.05, na.rm = TRUE),
  high_90 = quantile(VIP_Value, 0.95, na.rm = TRUE),
  CV_VIP = ifelse(mean(VIP_Value, na.rm = TRUE) == 0, NA, 
                  sd(VIP_Value, na.rm = TRUE) / mean(VIP_Value, na.rm = TRUE))
), by = .(Wavelength, Species)]

## subset species
# Display unique species names in the dataset
print(unique_species)

# Example: Select a subset of species names for plotting (modify this list as needed)
selected_species1 <- c("Acer", "Agonis", "Betula", "Claytosmunda", "Fagus", 
                      "Helianthus", "Myrica", "Osmunda", "Ostrya")
selected_species2 <- c("Phalaris", "Phragmites", "Populus", "Prunus", "Quercus", 
                       "Rubus", "Solidago", "Spiraea")

#Filter df_long to include only the selected species
vip_stats_subset <- vip_stats[Species %in% selected_species1]

# Plotting
ggplot(vip_stats_subset, aes(x = Wavelength)) +
  # 90% interval ribbon
  geom_ribbon(aes(ymin = low_90, ymax = high_90), fill = "blue", alpha = 0.5) +
  # 50% interval ribbon
  geom_ribbon(aes(ymin = low_50, ymax = high_50), fill = "blue", alpha = 0.9) +
  # Mean VIP line
  geom_line(aes(y = Mean_VIP, color = "Mean VIP"), linewidth = 0.6) +
  # CV VIP line
  geom_line(aes(y = CV_VIP, color = "CV VIP"), linewidth = 0.8) +
  labs(
    title = "Variable Importance in Projection: Herbarium Leaf PLSDA Model",
    x = "Wavelength (nm)",
    y = "VIP Values",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.5),
    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  # Scale x-axis to have ticks every 200 nm
  scale_x_continuous(breaks = seq(0, max(vip_stats$Wavelength, na.rm = TRUE), by = 200)) +
  scale_y_continuous() +
  # Custom colors for lines
  scale_color_manual(values = c("Mean VIP" = "black", "CV VIP" = "red")) +
  # Facet for each species
  facet_wrap(~ Species, scales = "free_y")

# Save the plot
ggsave("vip_plot_with_distributions.pdf", plot = p1, width = 8, height = 6)


#-------------------------------------------------------------------------------
#' @ LDA VIP Plots
#-------------------------------------------------------------------------------

# load data
vip_plsda <- fread(paste0(root_path,"/lda_varImp.csv"))

# Convert to data.table if needed
vip_lda <- as.data.table(vip_lda)

# Step 1: Rename "Variable" to "Wavelength" and remove backticks
setnames(vip_lda, "Variable", "Wavelength")
vip_lda[, Wavelength := as.numeric(gsub("`", "", Wavelength))]

# Step 2: Melt the data for long format
vip_lda_long <- melt(vip_lda, id.vars = c("Model", "Wavelength"), 
                     variable.name = "LD", 
                     value.name = "LD_Value")

# Step 3: Calculate statistics for each LDA component at each Wavelength
lda_stats <- vip_lda_long[, .(
  Mean_LDA = mean(LD_Value, na.rm = TRUE),
  low_50 = quantile(LD_Value, 0.25, na.rm = TRUE),
  high_50 = quantile(LD_Value, 0.75, na.rm = TRUE),
  low_90 = quantile(LD_Value, 0.05, na.rm = TRUE),
  high_90 = quantile(LD_Value, 0.95, na.rm = TRUE),
  CV_LDA = ifelse(mean(LD_Value, na.rm = TRUE) == 0, NA, 
                  sd(LD_Value, na.rm = TRUE) / mean(LD_Value, na.rm = TRUE))
), by = .(Wavelength, LD)]

# Step 4: Plot with facets for each LDA component
ggplot(lda_stats, aes(x = Wavelength)) +
  # 90% interval ribbon
  geom_ribbon(aes(ymin = low_90, ymax = high_90), fill = "blue", alpha = 0.5) +
  # 50% interval ribbon
  geom_ribbon(aes(ymin = low_50, ymax = high_50), fill = "blue", alpha = 0.9) +
  # Mean LDA line
  geom_line(aes(y = Mean_LDA, color = "Mean LDA"), linewidth = 0.8) +
  # CV LDA line
  geom_line(aes(y = CV_LDA, color = "CV LDA"), linewidth = 0.8) +
  labs(
    title = "Variable Importance in Projection: Herbarium Leaf LDA Model",
    x = "Wavelength (nm)",
    y = "LD VIP Values",
    color = "Metric"  # Legend title
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.5),
    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25),
    legend.position = "bottom",  # Legend at the bottom
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
  ) +
  # Scale x-axis to have ticks every 200 nm
  scale_x_continuous(breaks = seq(0, max(lda_stats$Wavelength, na.rm = TRUE), by = 200)) +
  scale_y_continuous() +
  # Manual color scale for lines
  scale_color_manual(values = c("Mean LDA" = "black", "CV LDA" = "red"))

# Save the plot
ggsave("vip_plot_with_distributions.pdf", plot = p1, width = 8, height = 6)

#-------------------------------------------------------------------------------
#' @Classification-Performance-Table
#-------------------------------------------------------------------------------

setwd("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/results/LDA/")

# Load data
perf <- readRDS("lda_performance_lda_training_sp10.rds")

# Convert to data frame
perf_df <- as.data.frame(perf$performance)

# Calculate mean and standard deviation for selected metrics
summary_stats <- perf_df %>%
  summarise(
    mean_accuracy = mean(Accuracy, na.rm = TRUE),
    sd_accuracy = sd(Accuracy, na.rm = TRUE),
    mean_precision = mean(Precision, na.rm = TRUE),
    sd_precision = sd(Precision, na.rm = TRUE),
    mean_balanced_accuracy = mean(`Balanced Accuracy`, na.rm = TRUE),
    sd_balanced_accuracy = sd(`Balanced Accuracy`, na.rm = TRUE)
  )

summary_stats <- summary_stats %>%
  mutate(
    Rank = "sp10",
    Split = "Training",
    Model = "LDA"
  )

lda_perf_testing_sp10 <- summary_stats
lda_perf_training_sp10 <- summary_stats

plsda_performance_plsda_training_species <- summary_stats
plsda_performance_plsda_training_genus <- summary_stats
plsda_performance_plsda_training_family <- summary_stats
plsda_performance_plsda_testing_species <- summary_stats
plsda_performance_plsda_testing_genus <- summary_stats
plsda_performance_plsda_testing_family <- summary_stats


# Combine LDA
combined_stats <- bind_rows(
  stats_lda_training_species,
  stats_lda_training_genus,
  stats_lda_testing_species,
  stats_lda_testing_genus,
  stats_lda_testing_family,
  stats_lda_training_family
)

# Combine PLS
combined_stats_plsda <- bind_rows(
  plsda_performance_plsda_training_species,
  plsda_performance_plsda_training_genus,
  plsda_performance_plsda_training_family,
  plsda_performance_plsda_testing_species,
  plsda_performance_plsda_testing_genus,
  plsda_performance_plsda_testing_family
)

# Combine sp10
combined_stats_sp10 <- bind_rows(
  lda_perf_testing_sp10,
  lda_perf_training_sp10
)


# Select and arrange columns in the specified order
final_table <- combined_stats_sp10 %>%
  select(Model, Rank, Split, mean_accuracy, sd_accuracy, mean_precision, sd_precision, mean_balanced_accuracy, sd_balanced_accuracy)

write.csv(final_table, "LDA_sp10_performance.csv")
