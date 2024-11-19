# Herbarium-Spectra  
Datasets and scripts for the **NPH paper**.

---

## Data 

Download Herbarium data here : https://drive.google.com/file/d/13kj_C4orE9QcR-01rbUHxUpX5ggX834J/view?usp=sharing

Download Kothari data here : https://drive.google.com/file/d/1MZBkpfeeSVQshd2gUXTuV-0XJbV3VWTb/view?usp=drive_link

These are normalized, 5nm resampled spectra

**`plotting/plotting_functions_spectra.R`** Plot means, quartiles, coefficient of variation of spectra.
  
---

## Data processing (completed by Dawson)
**`R/Spectrolab_process_data.R`**

Script for combining spectra and metadata files and manipulating spectra for output to analysis scripts.
  
---

## Trait Estimation  

The main script for trait estimation is:  

**`R/leaf_trait_prediction.R`**

This script takes as input either the **Kothari dataset** or the **herbarium dataset**, sets the metadata accordingly, and runs the auxiliary functions to estimate traits and generate output files for analysis.

#### Data Output  

The script generates the following files:  
- **`pls_*_split.rds`**  Testing and validation sample splits.  
- **`pls_*_segments.rds`**  Dataset iterations selecting random spectra per individual, based on the `accession` metadata column.  
- **`ncomp_onesigma-*.txt`**  Optimal number of components calculated using the one-sigma method from the `model_tune.R` function.  
- **`ncomp_otherMetrics-*.txt`**  Min or max values from PRESS, RMSEP, and R² statistics during cross-validation.  
- **`pls_*_opt_comp_models.csv`**  Cross-validation models from `model_tune`.  
- **`pls_*_final_models.rds`**  Final models across segments generated by the `model_build.R` function.  
- **`*performance.csv`**  Table containing statistics across accessions (nsamples, R², Bias, RMSE, %RMSE, Intercept, Slope)
- **`*obs-pred.csv`**  Data for plotting observed vs. predicted trait values and sample metadata (for plotting).
- **`*obs-pred-avg.csv`**  Similar to `*obs-pred.csv`, but averaged values are calculated for individuals.  
- **`pls_*_coefficients.csv`** Model coefficients across segments.  
- **`pls_*_vip.csv`** Variable importance in projection (VIP) scores across segments.  
- **`*plot.pdf`**  Plots showing ncomp optimization statistics.

### Auxiliary scripts
These are called by **`R/leaf_trait_prediction.R`**

### Plotting scripts
- **`R/plotting_functions_traits.R`** Observed vs predicted trait biplots, performance tables, VIP plots

---

## Taxonomic classification 

The main script for trait estimation is:  

**`R/leaf_classificationn.R`**

This script takes the **herbarium dataset** as input, sets the metadata accordingly, and runs the auxiliary functions to classify at different taxonomic levels using plsda and lda.

#### Data Output  
- **`classification_split.rds`**  Testing and validation sample splits.  
- **`classification_segments.rds`**  Dataset iterations selecting random spectra per individual, based on the `accession` metadata column.  
- **`classification_plsda_ncomp.txt`**  Optimal number of components calculated using the one-sigma method from the `model_tune.R` function.   
- **`classification_plsda_opt_comp_models.csv`**  Cross-validation models from `model_tune`.  
- **`classification_plsda_models.rds`**  Final models across segments generated by the `model_build.R` function.
- **`classification_lda_models.rds`**  Final models across segments generated by the `model_build.R` function.
- **`performance_*.csv`**  Table containing performance statistics across accessions.
- **`classification_*_coefficients.csv`** Model coefficients across segments.  
- **`classification_*_varImp.csv`** Variable importance in projection (VIP) scores across segments.
- **`CM_*.rds`** Data necessary for constructing confusion matrices using `plotting/plot_CM.R`.

### Auxiliary scripts
These are called by **`R/leaf_classification.R`**

### Plotting scripts
- **`plotting/plot_CM.R`** Confusion Matrices
- **`plotting/plotting_functions_classification.R`** VIP plots and performance tables.

---

## Analyses of herbarium factors on model success
- **`R/phylogenetic_diversity.R`**
Scripts for generating nearest taxon distance for regression models.
- **`regressions.R`**
Scripts for plotting classification results against herbarium factors, linear and logarthimic regressions, RF, and ANCOVA.
Still need to add scripts for analysis of trait predictions

