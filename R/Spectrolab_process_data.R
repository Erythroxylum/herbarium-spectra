################################################################################
# Load libraries
################################################################################
library("readr")
library("spectrolab")
library("shiny")

################################################################################
# Read data
################################################################################

# set wd
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/herbarium-spectra/")

# Spectral data within directory
data_refl = spectrolab::read_spectra(path = "~/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/Scan Data Archive/leaf_scans/",extract_metadata = TRUE)

####### get filenames
# Specify the path to the directory containing your spectra files
data_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/Scan Data Archive/leaf_scans/"

# List all the files in the directory
file_list <- list.files(path = data_dir, full.names = TRUE)

################################################################################
# Process spectra
################################################################################

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

####### vector normalize
data_refl_norm = spectrolab::normalize(data_refl)

#######
# Resample spectra with full width half maximum (FWHM) values

data_refl_resamp <- data_refl_norm

# Extract relevant data
values <- data_refl_resamp$value  # Spectral data
bands <- data_refl_resamp$bands  # Original wavelengths

# Define new wavelengths (350 to 2500 in 5 nm intervals)
new_wavelengths <- seq(350, 2500, by = 5)

# Initialize an empty list to store the resampled data
resampled_data <- list()

# Loop over each spectrum and resample
for (i in 1:nrow(values)) {
  # Extract the current spectrum
  sample_spectrum <- as.numeric(values[i, ])
  
  # Perform resampling with FWHM method
  resampled_sample <- prospectr::resample2(
    sample_spectrum,
    wav = as.numeric(bands),
    new.wav = new_wavelengths
  )
  
  # Add the resampled sample to the list
  resampled_data[[i]] <- resampled_sample
}

# Combine resampled data into a matrix
resampled_matrix <- do.call(rbind, resampled_data)

# Convert resampled data into the appropriate structure
#dimnames(resampled_matrix) <- list(NULL, as.character(new_wavelengths))

# Update the original data_refl object
data_refl_resamp$value <- resampled_matrix
data_refl_resamp$bands <- new_wavelengths

######### Trim ends
data_refl_resamp_norm_trim = data_refl_resamp[ , spectrolab::bands(data_refl_resamp, min = 400, max = 2400)]

# Trim ends,not normalized, plots
#data_refl_trim = data_refl[ , spectrolab::bands(data_refl, min = 400, max = 2400)]

################################################################################
# Match metadata and spectra
################################################################################

# set spectral data
spectra <- data_refl_resamp_norm_trim

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
#scans = sapply(split_names, `[`, 3)
#scan1 = strsplit(scans, split = ".", fixed = TRUE)
#scan        = sapply(scan1, `[`, 1)

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

# QC
unique(collector)
unique(class)
unique(order)
unique(family)
unique(genus)
unique(species)
unique(ddmmyyScanned)
hist(absoluteAge)
median(absoluteAge) # 92
plot(leafKg_m2)
plot(leafThickness)
unique(hasLMA) #  TRUE FALSE

unique(herbQuality) # "good"   "medium" "poor"
unique(damage) # "none"   "minor"  "medium" "major"
unique(glue) # FALSE TRUE
unique(leafStage) # "young"  "mature"



ggplot(database, aes(x = species, y = leafThickness)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Leaf Thickness by Species",
       x = "Species",
       y = "Leaf Thickness (mm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Troubleshoot:
# Identify which indices in species are NA
na_indices <- which(is.na(species))
na_indices
print(data_refl$name[na_indices])
print(data_refl$meta$time[na_indices])
print(accession_leaf[na_indices])

# Extract the relevant data from 'data_refl' based on 'na_indices'
bad_match_data <- data.frame(
  na_indices = na_indices,                            # First column: indices with NAs
  name = data_refl$name[na_indices],
  accession_leaf = accession_leaf[na_indices],
  time = data_refl$meta$time[na_indices]              # Third column: time values corresponding to those indices
)

# Write the data to a CSV file called "bad_match.csv"
write.csv(bad_match_data, file = "bad_match.csv", row.names = FALSE)




################################################################################
# Create table for analysis
################################################################################

full_data = data.frame(collector,accession, accession_leaf, leaf, scan, class, order, family, genus, species, growthForm,ddmmyyScanned,doy,absoluteAge,herbQuality,damage,glue,leafKg_m2,leafThickness,leafStage,as.matrix(spectra), check.names = F)

# Add green index
# Print the equation for the green index
cat("Green Index Equation: (ρ550 - ρ690) / (ρ550 + ρ690)\n")

# Calculate the green index for noResample
#greenIndex <- (full_data$`550.5` - full_data$`690.3`) / (full_data$`690.3` + full_data$`550.5`)

# Calculate the green index for 1nm resample
greenIndex <- (full_data$`550` - full_data$`690`) / (full_data$`690` + full_data$`550`)

# Add the greenIndex as the first column in the dataframe
full_data <- cbind(greenIndex, full_data)
# Optionally, rename the new column for clarity
colnames(full_data)[1] <- "greenIndex"

write.csv(full_data, file = "../dataHUH2024_sp25leaf563_norm_5nm_400-2400.csv")

# This does not insert a "filename" header in first row. Manually FIX and read in

unique(herbQuality) # "good"   "medium" "poor"
unique(damage) # "none"   "minor"  "medium" "major"
unique(glue) # FALSE TRUE
unique(leafStage) # "young"  "mature"







################################################################################
# Kothari data
################################################################################

# Load the required data
pressed <- readRDS("../../Kothari data/pressed_spec_unavg.rds")

# Normalize spectra
pressed_norm <- spectrolab::normalize(pressed)

# Assign data
names_col <- pressed_norm$names
meta_data <- pressed_norm$meta
bands <- pressed_norm$bands
values <- pressed_norm$value
values_df <- as.data.frame(values)
colnames(values_df) <- bands

## Resample to 5 nm across the full range (350–2500)
# Define the new wavelength interval
full_range <- seq(350, 2500, by = 5)

# Initialize an empty list to store the resampled data
resampled_data <- list()

# Loop over each row (sample) in the dataframe
for (i in 1:nrow(values_df)) {
  # Extract the current sample as a vector
  sample_spectrum <- as.numeric(values_df[i, ])
  
  # Perform resampling to 5 nm intervals
  resampled_sample <- prospectr::resample2(sample_spectrum, wav = as.numeric(colnames(values_df)), new.wav = full_range)
  
  # Add the resampled sample to the list
  resampled_data[[i]] <- resampled_sample
}

# Combine resampled data into a new data frame
resampled_df <- do.call(rbind, resampled_data)
colnames(resampled_df) <- full_range

# Trim resampled data to the range 400–2400
trimmed_wavelengths <- full_range[full_range >= 400 & full_range <= 2400]
resampled_df <- resampled_df[, as.character(trimmed_wavelengths)]

# Format species names to include only the first two words
meta_data$species <- sapply(strsplit(meta_data$Species, " "), function(x) paste(x[1:min(2, length(x))], collapse = " "))

# Remove rows for specific species
species_to_remove <- c("Parthenocissus quinquefolia", "Amelanchier humilis")
filtered_indices <- !meta_data$species %in% species_to_remove
resampled_df <- resampled_df[filtered_indices, ]
meta_data <- meta_data[filtered_indices, ]

# Add a new column to indicate TRUE/FALSE for PLSDA dataset
species_vector <- c("Quercus rubra", "Populus tremuloiddes", "Populus grandidentata", 
                    "Fagus grandifolia", "Betula populifolia", "Betula papyrifera", 
                    "Agonis flexuosa", "Acer saccharum", "Acer saccharinum", "Acer rubrum")
meta_data$sp10 <- meta_data$species %in% species_vector

# Add a new column for greenIndex
species_vector <- c("Quercus rubra", "Populus tremuloiddes", "Populus grandidentata", 
                    "Fagus grandifolia", "Betula populifolia", "Betula papyrifera", 
                    "Agonis flexuosa", "Acer saccharum", "Acer saccharinum", "Acer rubrum")
meta_data$greenIndex <- meta_data$species %in% species_vector

# Calculate the green index from wavelengths 550 and 690
greenIndex <- (resampled_df$`550` - resampled_df$`690`) / (resampled_df$`690` + resampled_df$`550`)

# Add the greenIndex as a new column to the combined data
combined_data <- cbind(name = names_col[filtered_indices], meta_data, greenIndex, resampled_df)

# rename ID column as 'accession'
colnames(combined_data)[colnames(combined_data) == "ID"] <- "accession"

# rename LMA
colnames(combined_data)[colnames(combined_data) == "LMA"] <- "leafKg_m2"

# Write the output to a CSV file
write.csv(combined_data, file = "../dataKothari_pressed_unavg_norm_5nm_400-2400.csv", row.names = FALSE)












################################################################################
# Create spectra object
################################################################################

#read in a data.frame
#specdf = read.table("train-test_clean_cocaSVC.txt", header=T, check.names = F)
revised_data <- read.csv("../fullDataHUH2024_sp25leaf636_noResample_400-2400.csv", header = T, check.names = F)
sp10 <- read.csv("../sp10HUH2024_sp25leaf636_noResample_400-2300.csv", header = T, check.names = F)
write.csv(sp10, file = "../sp10DataHUH2024_sp25leaf636_noResample_400-2300.csv")

#####
## implement quality filters
#####

data_herbQuality_noPoor <- revised_data_1[revised_data_1$herbQuality != "poor", ]
data_herbDam_noneMinor <- revised_data_1[revised_data_1$damage %in% c("none", "minor"), ]
data_noGlue <- revised_data_1[revised_data_1$glue %in% c("FALSE"), ]
data_leafMature <- revised_data_1[revised_data_1$leafStage %in% c("mature"), ]
cleanest_data <- revised_data_1[
  revised_data_1$herbQuality != "poor" &
    revised_data_1$damage %in% c("none", "minor") &
    revised_data_1$glue == "FALSE" &
    revised_data_1$leafStage == "mature",
]

# drop to 10 species:
sp10 <- revised_data[revised_data$species %in% c("Agonis flexuosa", "Acer rubrum", "Acer saccharinum", "Acer saccharum", "Betula papyrifera", "Betula populifolia", "Fagus grandifolia", "Populus grandidentata", "Populus tremuloides", "Quercus rubra"), ]


## all in 1
saveRDS(as_spectra(revised, name_idx = 1, meta_idxs = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,17,18)), file = "specHUH2024_1nm_sp25leaf639-1oct2024-herbQual_noPoor.rds")

saveRDS(as_spectra(data_herbDam_noneMinor, name_idx = 1, meta_idxs = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,17,18)), file = "specHUH2024_1nm_sp25leaf639-1oct2024-herbDam_noneMinor.rds")

saveRDS(as_spectra(data_noGlue, name_idx = 1, meta_idxs = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,17,18)), file = "specHUH2024_1nm_sp25leaf639-1oct2024-noGlue.rds")

saveRDS(as_spectra(data_leafMature, name_idx = 1, meta_idxs = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,17,18)), file = "specHUH2024_1nm_sp25leaf639-1oct2024-leafMature.rds")

saveRDS(as_spectra(cleanest_data, name_idx = 1, meta_idxs = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,17,18)), file = "specHUH2024_1nm_sp25leaf639-1oct2024-cleanest.rds")

saveRDS(as_spectra(sp10, name_idx = 1, meta_idxs = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,17,18)), file = "spectra/sp10HUH2024_sp25leaf636_noResample_400-2300.rds")


###########
### plot the number of leaves per species in each
###########

data_herbQuality_noPoor
data_herbDam_noneMinor
data_noGlue
data_leafMature
cleanest_data

# drop multiple scans per leaf
df_unique_full <- revised_data_1 %>%
  distinct(accession, .keep_all = TRUE)
df_unique_quality<- data_herbQuality_noPoor %>%
  distinct(accession, .keep_all = TRUE)
df_unique_damage <- data_herbDam_noneMinor %>%
  distinct(accession, .keep_all = TRUE)
df_unique_glue <- data_noGlue %>%
  distinct(accession, .keep_all = TRUE)
df_unique_mature <- data_leafMature %>%
  distinct(accession, .keep_all = TRUE)
df_unique_all <- cleanest_data %>%
  distinct(accession, .keep_all = TRUE)

# Add a new column to each dataframe to identify its source
df_unique_quality$source <- "Good or medium quality"
df_unique_damage$source <- "None or minor damage"
df_unique_glue$source <- "No Glue"
df_unique_mature$source <- "Mature leaves only"
df_unique_all$source <- "All filters"
df_unique_full$source <- "No filters"

# Combine all species columns into one dataframe
combined_data <- rbind(
  df_unique_full[, c("species", "source")],
  df_unique_mature[, c("species", "source")],
  df_unique_glue[, c("species", "source")],
  df_unique_quality[, c("species", "source")],
  df_unique_damage[, c("species", "source")],
  df_unique_all[, c("species", "source")]
)

# Count the number of rows per species and source
library(dplyr)

# Build the species counts table, maintaining the order of combined_data$source
species_counts <- combined_data %>%
  count(source, species) %>%
  mutate(source = factor(source, levels = unique(combined_data$source))) %>%
  arrange(source)

# Rename the columns for clarity
colnames(species_counts) <- c("source", "species", "count")

# Calculate the total number of rows per source for labeling the top of the bars
total_counts <- aggregate(count ~ source, data = species_counts, sum)

# Plot the stacked bar chart
library(viridisLite)
ggplot(species_counts, aes(x = source, y = count, fill = species)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option = "turbo") + # Use the turbo palette from viridisLite
  labs(x = "Quality Filters", y = "Number of Individuals", title = "Number of Individuals Scanned Across Datasets") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis text by 45 degrees












####### old

## The sample names are in column 1. Mark other metadata columns
spec_1 = as_spectra(revised_data_1, name_idx = 1, meta_idxs = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,17,18))
spec_5 = as_spectra(revised_data_5, name_idx = 1, meta_idxs = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,17,18))
## save the spectral object as RDS to be loaded in plsda software
saveRDS(spec_1, file = "specHUH2024_1nm_sp25leaf639-1oct2024-full.rds")
saveRDS(spec_5, file = "specHUH2024_5nm_sp25leaf639-1oct2024-full.rds")

## And now you have a spectra object with sample names and metadata...
spec
spec["81-32-854",]
spec["spec_12049"]



