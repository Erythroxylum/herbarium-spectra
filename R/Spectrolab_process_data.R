#'------------------------------------------------------------------------------
#' @title Spectral raw data processing
#'------------------------------------------------------------------------------

#' @description Scripts to input spectra files or spectra objects and metadata
#' and export matrixes in csv. This is customized for the HUH herbarium
#' dataset and the Kothari pressed-leaf dataset.
#' Plotting functions for datasets can be found at the bottom.
#' 
#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------
library(readr)
library(spectrolab)
library(shiny)
library(dplyr)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# set wd
setwd("~/Downloads/")

#'------------------------------------------------------------------------------
#' @Process-Herbarium-Spectra
#-------------------------------------------------------------------------------

# Spectral data within directory
data_refl = spectrolab::read_spectra(path = "~/Downloads/DMWhiteHUHspec1_raw_spectra_files/",extract_metadata = TRUE)

####### get filenames
# Specify the path to the directory containing your spectra files
data_dir <- "~/Downloads/DMWhiteHUHspec1_raw_spectra_files/"

# List all the files in the directory
file_list <- list.files(path = data_dir, full.names = TRUE)

###############################
#  match sensor overlap

## Guess "good" splicing bands
#splice_guess = spectrolab::guess_splice_at(data_refl)
## Viz the spectra and confirm that the splicing bands make sense
#spectrolab::plot_interactive(data_refl)

# Pick the final splicing bands. Could be the same as `splice_guess`
# Dawson did this by selecting the right side of the break. Worked.
splice_at = c(991.3, 1902.5)

# Finally, match the sensor overlap, interpolate_wv=extent around splice_at values over which the splicing factors will be calculated.
data_refl = spectrolab::match_sensors(x = data_refl, splice_at = splice_at, interpolate_wvl = c(5, 1))

###############################
# Normalize spectra

data_refl_norm = spectrolab::normalize(data_refl)

###############################
# Match metadata and spectra

# set spectral data: data_refl or data_refl_norm
spectra <- data_refl
#spectra <- data_refl_norm

# Extract just the base filenames (without directory path and extension)
names(spectra) <- tools::file_path_sans_ext(basename(file_list))

# Read metadata spreadsheet
metadata <- readr::read_csv("~/Downloads/DMWhiteHUHspec1_sp25leaf560_metadata.csv")

# Get spectra names (without file extension)
spectra_names <- names(spectra)
spectra_base <- tools::file_path_sans_ext(basename(spectra_names))

# Match metadata rows using filename (no extension)
metadata_base <- tools::file_path_sans_ext(metadata$filename)
match_idx <- match(spectra_base, metadata_base)

# Reorder metadata to match spectra
metadata_matched <- metadata[match_idx, ]

# Combine metadata and spectra into one data frame
spectra_matrix <- as.data.frame(as.matrix(spectra))  # convert spectra to numeric matrix
full_data <- cbind(metadata_matched, spectra_matrix)

# Write to CSV
write.csv(full_data, file = "../DMWhiteHUHspec1_sp25leaf560_noResample_350-2500.csv")


#-------------------------------------------------------------------------------
#' @Plot-spectra
#-------------------------------------------------------------------------------

# Read spectra
root_path <- getwd()
frame <- fread(paste0(root_path, "/DMWhiteHUHspec1_sp25leaf560_noResample_350-2500.csv"))

# Define bands of interest
#bands <- seq(450, 2400, by = 5)
# Identify column index where spectral data starts
spectral_start <- which(colnames(frame) == "338.6")

# Subset metadata columns explicitly
meta <- frame[, c("specimenIdentifier", "scientificName", "Genus", "Family", "growthForm")]

# Subset spectral columns from 338.6 to the end
spectra <- frame[, spectral_start:ncol(frame)]
cbands <- as.character(bands)

## Subset data
# Define species of interest
species_of_interest <- "Populus tremuloides"

# Filter metadata and spectra by species
filtered_sp_meta <- meta[meta$scientificName == species_of_interest, ]
filtered_sp_spectra <- spectra[meta$scientificName == species_of_interest, ]

# Add a sample identifier column
filtered_sp_meta$sample_id <- seq_len(nrow(filtered_sp_meta))

# Combine metadata and spectra
combined_data <- cbind(filtered_sp_meta, filtered_sp_spectra)

# Identify spectral column names dynamically
spectral_cols <- colnames(filtered_sp_spectra)

# Reshape to long format for plotting
library(tidyr)
long_data <- combined_data %>%
  pivot_longer(
    cols = all_of(spectral_cols),
    names_to = "wavelength",
    values_to = "reflectance"
  )

# Convert wavelength to numeric
long_data$wavelength <- as.numeric(long_data$wavelength)

# Plot
library(ggplot2)
ggplot(long_data, aes(x = wavelength, y = reflectance, color = specimenIdentifier, group = sample_id)) +
  geom_line(alpha = 0.7) +
  labs(
    title = paste("Spectra for", species_of_interest),
    x = "Wavelength (nm)",
    y = "Reflectance",
    color = "Specimen"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 10, face = "bold")
  )

# Save the plot
ggsave("spectra_plot_Poptre_all.pdf", plot = last_plot(), device = "pdf", width = 5, height = 3.5)


#-------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------

