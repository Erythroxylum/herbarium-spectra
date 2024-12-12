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
library("readr")
library("spectrolab")
library("shiny")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

# set wd
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/herbarium-spectra/")

#'------------------------------------------------------------------------------
#' @Read-Data
#-------------------------------------------------------------------------------

# Spectral data within directory
data_refl = spectrolab::read_spectra(path = "~/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/Scan Data Archive/leaf_scans/",extract_metadata = TRUE)

####### get filenames
# Specify the path to the directory containing your spectra files
data_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/Scan Data Archive/leaf_scans/"

# List all the files in the directory
file_list <- list.files(path = data_dir, full.names = TRUE)

#'------------------------------------------------------------------------------
#' @Process-Spectra-HUH
#-------------------------------------------------------------------------------
#######
# match sensor overlap

## Guess "good" splicing bands
#splice_guess = spectrolab::guess_splice_at(data_refl)
## Viz the spectra and confirm that the splicing bands make sense
#spectrolab::plot_interactive(data_refl)

# Pick the final splicing bands. Could be the same as `splice_guess`
# Dawson did this by selecting the right side of the break. Worked.
splice_at = c(991.3, 1902.5)

# Finally, match the sensor overlap, interpolate_wv=extent around splice_at values over which the splicing factors will be calculated.
data_refl = spectrolab::match_sensors(x = data_refl, splice_at = splice_at, interpolate_wvl = c(5, 1))

#-------------------------------------------------------------------------------
#' @Normalize
#-------------------------------------------------------------------------------

# Normalize spectra
data_refl_norm = spectrolab::normalize(data_refl)

#-------------------------------------------------------------------------------
#' @Match-metadata-and-spectra
#-------------------------------------------------------------------------------

# set spectral data, normalized or otherwise
spectra <- data_refl_norm

# rename spectra
#names(data_refl) = meta(x = data_refl, label = "name", simplify = TRUE)
###### redo to get filenames
# Extract just the base filenames (without directory path and extension)
names(spectra) <- tools::file_path_sans_ext(basename(file_list))
######

# Parse names of spectra files to get accession number information
split_names = strsplit(names(spectra), split = "_", fixed = F) # split spectra filenames to get accession
accession   = sapply(split_names, `[`, 1)
leaf        = sapply(split_names, `[`, 2)
scan <- sapply(strsplit(sapply(split_names, `[`, 3), split = ".", fixed = TRUE), `[`, 1)
accession_leaf = paste(accession,"_",leaf,sep="")

# read metadata file
database <- readr::read_csv("/Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/dataSheet_LMA_metadata-sp25leaf639-1oct2024.csv")

# get accesssion_leaf number to match rows of database to data_refl
match_idx = match(accession_leaf, database$catalogNumber_leaf)

# Pull other columns from metadata file:
collector <- database$Collector[match_idx]
class <- database$Class[match_idx]
order <- database$Order[match_idx]
family <- database$Family[match_idx]
genus <- database$Genus[match_idx]
species <- database$species[match_idx]
ddmmyyScanned <- database$ddmmyyScanned[match_idx]
absoluteAge <- database$absoluteAge[match_idx]
herbQuality <- database$herbQuality[match_idx]
damage <- database$damage[match_idx]
glue <- database$glue[match_idx]
leafKg_m2 <- database$leafKg_m2[match_idx]
leafThickness <- database$leafThickness[match_idx]
leafStage <- database$leafStage[match_idx]
hasLMA <- database$hasLMA[match_idx]
doy <- database$doy[match_idx]
growthForm <- database$growthForm[match_idx]

#-------------------------------------------------------------------------------
#' @Build-and-write-data-table
#-------------------------------------------------------------------------------

full_data = data.frame(collector,accession, accession_leaf, leaf, scan, class, order, family, genus, species, growthForm,ddmmyyScanned,doy,absoluteAge,herbQuality,damage,glue,leafKg_m2,leafThickness,leafStage,as.matrix(spectra), check.names = F)

## Add green index

# Calculate the green index for noResample
greenIndex <- (full_data$`550.5` - full_data$`690.3`) / (full_data$`690.3` + full_data$`550.5`)

# Add the greenIndex as the first column in the dataframe
full_data <- cbind(greenIndex, full_data)

# Optionally, rename the new column for clarity
colnames(full_data)[1] <- "greenIndex"

write.csv(full_data, file = "../dataHUH2024_sp25leaf563_norm_350-2500.csv")

# This does not insert a "filename" header in first row. Manually FIX and read in









#'------------------------------------------------------------------------------
#' @Kothari-Data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#' @Load-Spectra-Object
#-------------------------------------------------------------------------------

pressed <- readRDS("../../Kothari data/pressed_spec_unavg_v3.rds")

#-------------------------------------------------------------------------------
#' @Splice-overlap-region
#-------------------------------------------------------------------------------

## Guess "good" splicing bands
#splice_guess = spectrolab::guess_splice_at(data_refl)
## Viz the spectra and confirm that the splicing bands make sense
#spectrolab::plot_interactive(pressed)

## Not doing splicing for these analyses, only removing splice regions.
# Pick the final splicing bands. Could be the same as `splice_guess`
#splice_at = c(995)
# Finally, match the sensor overlap, interpolate_wv=extent around splice_at values over which the splicing factors will be calculated.
#data_refl = spectrolab::match_sensors(x = pressed, splice_at = splice_at, interpolate_wvl = c(5, 1))

#-------------------------------------------------------------------------------
#' @Normalize-optional
#-------------------------------------------------------------------------------

# Normalize spectra
pressed_norm <- spectrolab::normalize(pressed)

# Assign data
names_col <- pressed$names
meta_data <- pressed$meta
bands <- pressed$bands
values <- pressed$value
values_df <- as.data.frame(values)
colnames(values_df) <- bands

#-------------------------------------------------------------------------------
#' @Revise-metadata
#-------------------------------------------------------------------------------

# Format species names to include only the first two words
meta_data$species <- sapply(strsplit(meta_data$Species, " "), function(x) paste(x[1:min(2, length(x))], collapse = " "))

# Remove rows for specific species
species_to_remove <- c("Parthenocissus quinquefolia", "Amelanchier humilis")
filtered_indices <- !meta_data$species %in% species_to_remove
values_df <- values_df[filtered_indices, ]
meta_data <- meta_data[filtered_indices, ]

# Add a new column to indicate TRUE/FALSE for PLSDA dataset
species_vector <- c("Quercus rubra", "Populus tremuloiddes", "Populus grandidentata", 
                    "Fagus grandifolia", "Betula populifolia", "Betula papyrifera", 
                    "Agonis flexuosa", "Acer saccharum", "Acer saccharinum", "Acer rubrum")
meta_data$sp10 <- meta_data$species %in% species_vector

# Calculate the green index from wavelengths 550 and 690
greenIndex <- (values_df$`550` - values_df$`690`) / (values_df$`690` + values_df$`550`)

# Add the greenIndex as a new column to the combined data
combined_data <- cbind(name = names_col[filtered_indices], meta_data, greenIndex, values_df)

# rename ID column as 'accession'
colnames(combined_data)[colnames(combined_data) == "ID"] <- "accession"

# rename LMA
colnames(combined_data)[colnames(combined_data) == "LMA"] <- "leafKg_m2"

# rename growthForm
colnames(combined_data)[colnames(combined_data) == "GrowthForm"] <- "growthForm"

# Write the output to a CSV file
write.csv(combined_data, file = "../dataKothari_pressed_unavg_noResample_350-2500.csv", row.names = FALSE)








#-------------------------------------------------------------------------------
#' @Plot-spectra
#-------------------------------------------------------------------------------

# Read spectra
root_path <- getwd()
frame <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf563_ref5nm_norm_450-2400.csv"))

# Define bands of interest
bands <- seq(450, 2400, by = 5)
cbands <- as.character(bands)

# Extract meta and spectra
meta <- frame[, c("accession", "species", "genus", "family", "growthForm")]
spectra <- frame[, ..cbands]

## Subset data
# Define species of interest
species_of_interest <- "Populus tremuloides"
# Filter meta
filtered_sp_meta <- meta[species == species_of_interest]
# Filter spectra
filtered_sp_spectra <- spectra[meta$species == species_of_interest]

# Select Accessions of Interest
# E.g.pick the first 3 unique accessions
accessions_of_interest <- unique(filtered_sp_meta$accession)[20:29]
#or make a vector of accession numbers
#accessions_of_interest <- c("")
# Filter meta and spectra for the selected accessions
filtered_meta <- filtered_sp_meta[accession %in% accessions_of_interest]
filtered_spectra <- filtered_sp_spectra[filtered_sp_meta$accession %in% accessions_of_interest]

# Add a sample identifier to differentiate lines
filtered_meta$sample_id <- seq_len(nrow(filtered_meta))

# Combine filtered_meta and filtered_spectra
combined_data <- cbind(filtered_meta, filtered_spectra)

# Reshape the combined data to long format for plotting
long_data <- combined_data %>%
  tidyr::pivot_longer(
    cols = all_of(cbands), # Use column names corresponding to wavelengths
    names_to = "wavelength",                 # New column for wavelengths
    values_to = "reflectance"                # New column for reflectance values
  )

# Convert wavelength to numeric for correct plotting
long_data$wavelength <- as.numeric(long_data$wavelength)

# Plot the spectra using ggplot
ggplot(long_data, aes(x = wavelength, y = reflectance, color = accession, group = sample_id)) +
  geom_line(alpha = 0.7) +
  labs(
    title = paste("Spectra for", species_of_interest, ":724305 problem"),
    x = "Wavelength (nm)",
    y = "Reflectance",
    color = "Accession"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 10, face = "bold")
  )

ggsave("spectra_plot_Poptre_problem724305_sm.pdf", plot = last_plot(), device = "pdf", width = 5, height = 3.5)


#-------------------------------------------------------------------------------
#' @End
#-------------------------------------------------------------------------------

