# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# do cwtnorm then ref
# do each HUH trait for each ggplot and then make final plot together.

# load Pressed data
Kothari_data <- read.csv("data/dataKothari_pressed_unavg_cwt5nm_norm_450-2400.csv")

# set trait:  nothing or leafKg_m2 for LMA, C, Ca, car, cel or cellulose, chlA, N, sol or solubles
trait_name <- "solubles"

# Step 1: Load HUH trait data
HUH_data <- read.csv(paste0("results/model_transfer/Coef-Kothari_Spec-HUH_refnorm", trait_name, "_1350-2400/Coef-Kothari_Spec-HUH_refnorm", trait_name, "_1350-2400_predicted_trait.csv"))

# Step 2: parse HUH data
HUH_parse <- HUH_data %>%
  select(species, predicted_trait_mean) %>%
  mutate(dataset = "Herbarium", trait_value = predicted_trait_mean) %>%
  select(species, trait_value, dataset)

# Step 3: parse species in HUH_data
HUH_species <- HUH_parse %>% select(species) %>% distinct()

# Step 4: parse Kothari trait data
Kothari_filtered <- Kothari_data %>%
  filter(species %in% HUH_species$species) %>%
  select(species, trait_value = all_of(trait_name)) %>% # Dynamically select the trait column
  mutate(dataset = "Pressed") %>% 
  select(species, trait_value, dataset)

# Step 3: Combine the two datasets
combined_data_refnorm_LMA <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_refnorm_C <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_refnorm_Ca <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_refnorm_car <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_refnorm_cel <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_refnorm_chlA <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_refnorm_N <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_refnorm_sol <- bind_rows(HUH_parse, Kothari_filtered)


# Step 3: Combine the two datasets
#combined_data_cwtnorm_LMA <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_cwtnorm_C <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_cwtnorm_Ca <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_cwtnorm_car <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_cwtnorm_cel <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_cwtnorm_chlA <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_cwtnorm_N <- bind_rows(HUH_parse, Kothari_filtered)
combined_data_cwtnorm_sol <- bind_rows(HUH_parse, Kothari_filtered)



## Step 4: Create a split violin plot for trait value distributions
# set combined data and scale_y_continuous name

# plots with no x axis labels
trait_LMA <- ggplot(combined_data_refnorm_LMA, aes(x = species, y = trait_value, fill = dataset)) +
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

trait_C <- ggplot(combined_data_refnorm_C, aes(x = species, y = trait_value, fill = dataset)) +
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

trait_Ca <- ggplot(combined_data_refnorm_Ca, aes(x = species, y = trait_value, fill = dataset)) +
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

trait_car <- ggplot(combined_data_refnorm_car, aes(x = species, y = trait_value, fill = dataset)) +
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

trait_cel <- ggplot(combined_data_refnorm_cel, aes(x = species, y = trait_value, fill = dataset)) +
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

trait_chlA <- ggplot(combined_data_refnorm_chlA, aes(x = species, y = trait_value, fill = dataset)) +
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

trait_N <- ggplot(combined_data_refnorm_N, aes(x = species, y = trait_value, fill = dataset)) +
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
trait_sol <- ggplot(combined_data_refnorm_sol, aes(x = species, y = trait_value, fill = dataset)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE, scale = "width") + # Adjust scale for wider violins
  stat_summary(fun.data = "mean_se", geom = "pointrange", 
               position = position_dodge(0.8), show.legend = FALSE, size = 0.1) + # Adjust size for visibility
  scale_x_discrete(labels = function(x) lapply(x, function(label) bquote(italic(.(label))))) + # Italicize species names
  scale_y_continuous(name = expression(Solubles ~ (mg ~ g^{-1}))) + # Correct y-axis for trait_value (continuous)
  scale_fill_manual(name = "Dataset: Normalized,1350-2400 nm, ", values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
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
library(patchwork)
combined_plot <- trait_LMA / trait_C / trait_Ca / trait_car / trait_cel / trait_chlA / trait_N / trait_sol

ggsave("results/Figure-trait_violins_refnorm1350.pdf", plot = combined_plot, width = 8.5, height = 11)








###############################
# Violin
Vplot <- ggplot(combined_data, aes(x = trait_value, y = species, fill = dataset)) +
  geom_violin(scale = "width", trim = TRUE) +
  labs(
    title = "Trait Value Distributions by Species: C",
    x = "Trait Value",
    y = "Species",
    fill = "Dataset"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Herbarium" = "#56B4E9", "Pressed" = "#E69F00")) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "bottom"
  )

## boxplot
ggplot(combined_data, aes(x = species, y = trait_value, fill = dataset)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.7, outlier.size = 0.8) +
  scale_x_discrete(name = "Species") +
  scale_y_continuous(name = "Trait Value") +
  scale_fill_manual(values = c("Herbarium (Predicted)" = "#56B4E9", "Pressed (Measured)" = "#E69F00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = "Trait Value Distributions by Species",
    fill = "Dataset"
  )



