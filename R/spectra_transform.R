#'------------------------------------------------------------------------------
#' @title Spectral homogenization among sensors
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Library
#-------------------------------------------------------------------------------

library(data.table)
library(CWT)
library(spectrolab)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- getwd()

#'------------------------------------------------------------------------------
#' @HUH
#-------------------------------------------------------------------------------

# Read database
data <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf561_noResample_350-2500.csv"))
data <- fread(paste0(root_path, "/data/dataHUH2024_sp25leaf561_norm_350-2500.csv"))

# Get meta and spectra
meta <- data[, .SD, .SDcols = 1:22]
spectra <- data[, .SD, .SDcols = 23:ncol(data)]

# Difference between consecutive bands (resolution)
bands <- as.numeric(colnames(spectra))
resolution <- diff(bands)
range(resolution) # Lower and higher resolution
res <- 5.0 # Define resolution to use

# Resample spectra
resampled_bands <- ceiling(seq(340, 2510, res))
FHWM <- rep(res, length(resampled_bands))

spectra_resampled <- resampling_FWHM(spectra = spectra,
                                     wavelengths = bands,
                                     new_wavelengths = resampled_bands,
                                     FWHM = FHWM)

plot(as.numeric(colnames(spectra_resampled)), 
     as.matrix(spectra_resampled[50,])[1,])

# Remove edges
spectra_resampled <- spectra_resampled[, .SD, .SDcols = 3:(ncol(spectra_resampled)-2)]

plot(as.numeric(colnames(spectra_resampled)), 
     as.matrix(spectra_resampled[50,])[1,])

## Trim edges <450, >2400
spectra_ref_export <- spectra_resampled[, .SD, .SDcols = 21:411]

fwrite(cbind(meta, spectra_ref_export), 
       paste0(root_path, "/data/dataHUH2024_sp25leaf561_ref5nm_norm_450-2400.csv"))

# Transform spectra
spectra_transformed <- cwt(t = spectra_resampled, 
                           scales = c(2, 3, 4), 
                           variance = 1, 
                           summed_wavelet = TRUE, 
                           threads = 1L)

plot(as.numeric(colnames(spectra_transformed)), 
     as.matrix(spectra_transformed[50,])[1,])
abline(v = 450)
abline(v = 2400)

# Remove edges <450, >2400
spectra_cwt_export <- spectra_transformed[, .SD, .SDcols = 21:411]

# Export
fwrite(cbind(meta, spectra_cwt_export), 
       paste0(root_path, "/data/dataHUH2024_sp25leaf561_cwt5nm_norm_450-2400.csv"))



#'------------------------------------------------------------------------------
#' @Kothari
#-------------------------------------------------------------------------------

# Read database
data <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_noResample_350-2500.csv"))
data <- fread(paste0(root_path, "/data/dataKothari_pressed_unavg_noResample_norm_350-2500.csv"))

# Get meta and spectra
meta <- data[, .SD, .SDcols = 1:60]
spectra <- data[, .SD, .SDcols = 61:ncol(data)]

# Difference between consecutive bands (resolution)
bands <- as.numeric(colnames(spectra))
resolution <- diff(bands)
range(resolution) # Lower and higher resolution
res <- 5 # Define resolution to use

# Resample spectra
resampled_bands <- ceiling(seq(350, 2500, res))
FHWM <- rep(res, length(resampled_bands))

spectra_resampled <- resampling_FWHM(spectra = spectra,
                                     wavelengths = bands,
                                     new_wavelengths = resampled_bands,
                                     FWHM = FHWM)

plot(as.numeric(colnames(spectra_resampled)), 
     as.matrix(spectra_resampled[50,])[1,])


# Remove edges
spectra_resampled <- spectra_resampled[, .SD, .SDcols = 3:(ncol(spectra_resampled)-2)]

plot(as.numeric(colnames(spectra_resampled)), 
     as.matrix(spectra_resampled[6,])[1,])
abline(v=450)
abline(v=2300)

## Trim edges <450, >2400
kothari_ref_export <- spectra_resampled[, .SD, .SDcols = 19:409]

fwrite(cbind(meta, kothari_ref_export), 
       paste0(root_path, "/data/dataKothari_pressed_unavg_ref5nm_norm_450-2400.csv"))

# Transform spectra
spectra_transformed <- cwt(t = spectra_resampled, 
                           scales = c(2, 3, 4), 
                           variance = 1, 
                           summed_wavelet = TRUE, 
                           threads = 1L)

plot(as.numeric(colnames(spectra_transformed)), 
     as.matrix(spectra_transformed[50,])[1,])
abline(v = 450)
abline(v = 2400)

## Trim edges <450, >2400
kothari_cwt_export <- spectra_transformed[, .SD, .SDcols = 19:409]

# Export
fwrite(cbind(meta, kothari_cwt_export), 
       paste0(root_path, "/data/dataKothari_pressed_unavg_cwt5nm_450-2400.csv"))
