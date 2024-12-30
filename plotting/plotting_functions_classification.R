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

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#' @Classification-Performance-Table-S3
#-------------------------------------------------------------------------------

library(dplyr)

# Define the file paths and metadata
files <- list(
  list(path = "out_classify/performance_lda_testing_sp10.rds", dataset = "Herbarium", rank = "species", model = "LDA", samples = "10spp"),
  list(path = "out_classify/performance_plsda_testing_sp10.rds", dataset = "Herbarium", rank = "species", model = "PLSDA", samples = "10spp"),
  list(path = "out_classify_genus/performance_plsda_testing_sp10.rds", dataset = "Herbarium", rank = "genus", model = "PLSDA", samples = "6genera"),
  list(path = "out_classify_genus/performance_lda_testing_sp10.rds", dataset = "Herbarium", rank = "genus", model = "LDA", samples = "6genera"),
  list(path = "out_classify/performance_lda_testing.rds", dataset = "Herbarium", rank = "species", model = "LDA", samples = "25spp"),
  list(path = "out_classify/performance_plsda_testing.rds", dataset = "Herbarium", rank = "species", model = "PLSDA", samples = "25spp"),
  list(path = "out_classify_genus/performance_plsda_testing.rds", dataset = "Herbarium", rank = "genus", model = "PLSDA", samples = "17genera"),
  list(path = "out_classify_genus/performance_lda_testing.rds", dataset = "Herbarium", rank = "genus", model = "LDA", samples = "17genera"),
  list(path = "out_classify_kot/performance_plsda_testing_sp10.rds", dataset = "Pressed", rank = "species", model = "PLSDA", samples = "10spp"),
  list(path = "out_classify_kot/performance_lda_testing_sp10.rds", dataset = "Pressed", rank = "species", model = "LDA", samples = "10spp"),
  list(path = "out_classify_kot_LatinGenus/performance_plsda_testing_sp10.rds", dataset = "Pressed", rank = "genus", model = "PLSDA", samples = "6genera"),
  list(path = "out_classify_kot_LatinGenus/performance_lda_testing_sp10.rds", dataset = "Pressed", rank = "genus", model = "LDA", samples = "6genera")
)

# Initialize an empty data frame for combined statistics
combined_stats <- data.frame()

# Loop through the files and process each
for (file_info in files) {
  # Load the performance data
  perf <- readRDS(file_info$path)
  
  # Convert to data frame
  perf_df <- as.data.frame(perf$performance)
  
  # Calculate summary statistics
  summary_stats <- perf_df %>%
    summarise(
      mean_accuracy = mean(Accuracy, na.rm = TRUE),
      sd_accuracy = sd(Accuracy, na.rm = TRUE),
      mean_precision = mean(Precision, na.rm = TRUE),
      sd_precision = sd(Precision, na.rm = TRUE),
      mean_balanced_accuracy = mean(`Balanced Accuracy`, na.rm = TRUE),
      sd_balanced_accuracy = sd(`Balanced Accuracy`, na.rm = TRUE)
    ) %>%
    mutate(
      Dataset = file_info$dataset,
      Rank = file_info$rank,
      Model = file_info$model,
      Samples = file_info$samples,
      ncomp = ncomp
    )
  
  # Append to the combined stats table
  combined_stats <- bind_rows(combined_stats, summary_stats)
}

# Select and arrange columns in the specified order
final_table <- combined_stats %>%
  select(Rank, Dataset, Model, Samples, mean_accuracy, sd_accuracy, mean_precision, sd_precision, mean_balanced_accuracy, sd_balanced_accuracy)

# Save the final table to a CSV file
write.csv(final_table, "Table_classification_performance.csv", row.names = FALSE)






#-------------------------------------------------------------------------------
#' @PLSDA_VIP_Plots_by_species
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
ggsave("vip_plot_herb_plsda_25spp.pdf", plot = vipsp, width = 8.5, height = 9.5)



#-------------------------------------------------------------------------------
#' @PLSDA_VIP_Plots_across_datasets
#-------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(gridExtra)

# Define the file paths and metadata
file_info <- list(
  list(path = "out_classify/varImp_plsda_sp10.csv", dataset = "Herbarium", rank = "species", model = "PLSDA", samples = "10 spp"),
  list(path = "out_classify_kot/varImp_plsda_sp10.csv", dataset = "Pressed", rank = "species", model = "PLSDA", samples = "10 spp"),
  list(path = "out_classify_genus/varImp_plsda_sp10.csv", dataset = "Herbarium", rank = "genus", model = "PLSDA", samples = "6 genera"),
  list(path = "out_classify_kot_LatinGenus/varImp_plsda_sp10.csv", dataset = "Pressed", rank = "genus", model = "PLSDA", samples = "6 genera"),
  list(path = "out_classify/varImp_plsda.csv", dataset = "Herbarium", rank = "species", model = "PLSDA", samples = "25 spp"),
  list(path = "out_classify_genus/varImp_plsda.csv", dataset = "Herbarium", rank = "genus", model = "PLSDA", samples = "17 genera")
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
    scale_x_continuous(breaks = seq(0, max(vip_stats$Wavelength, na.rm = TRUE), by = 200))
  
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
pdf("Fig-S_vip_classification_pressed_vs_herb.pdf", width = 7, height = 9)

# First page with six plots
grid.arrange(grobs = plots[1:6], ncol = 1)
dev.off()

# Second page with two plots, 8 x 5 layout
#pdf("vip_plots2.pdf", width = 8, height = 4)
#grid.arrange(grobs = plots[5:6], ncol = 1)
#dev.off()


