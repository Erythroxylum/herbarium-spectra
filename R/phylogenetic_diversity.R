#'------------------------------------------------------------------------------
#' @title Get Phylogenetic Distances
#'------------------------------------------------------------------------------

#' @description Estimates distances from Smith & Brown tree using the 
#' V.PhyloMaker2 package
#' @return Return distance matrices as .rds files.

################################################################################
#' @warming If you have to many communities this need to be run using HPC
################################################################################

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(phytools)
library(V.PhyloMaker2)
#library(picante)
library(bit64)

#'------------------------------------------------------------------------------
#' Source code
setwd("/panfs/jay/groups/17/cavender/guzman/MMDB")
#source("R/beta_rcpp.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param SPCD
#' @param out_path
#' @param threads
#' @param progress

#'------------------------------------------------------------------------------
#' Function

#phylogenetic_diversity <- function(csv_path, threads = 1) {
  
  # Read the dataset from a CSV file
  data <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Spectroscopy/HUH-Spectra-2024/RAnalysesHUHSpectra/fullDataHUH2024_sp25leaf636_noResample_400-2300.csv")
  
  # Extract relevant columns and retain unique species
  species_data <- unique(data[, c("species", "genus", "family")])
  #species_data$species <- gsub(" ", "_", species_data$species)
  
  species_data <- data.frame(
    species = c("Myrica gale", "Acer rubrum", "Betula papyrifera", "Betula populifolia",
                "Solidago gigantea", "Populus tremuloides", "Populus grandidentata", 
                "Prunus pensylvanica", "Osmunda claytoniana", "Phragmites australis",
                "Phalaris arundinacea", "Acer saccharinum", "Rubus odoratus", "Spiraea latifolia", 
                "Ostrya virginiana", "Acer saccharum", "Acer spicatum", "Quercus rubra", 
                "Rubus idaeus", "Osmunda regalis", "Fagus grandifolia", "Prunus serotina", 
                "Helianthus divaricatus", "Solidago altissima", "Agonis flexuosa"),
    genus = c("Myrica", "Acer", "Betula", "Betula", "Solidago", "Populus", "Populus", 
              "Prunus", "Osmunda", "Phragmites", "Phalaris", "Acer", "Rubus", "Spiraea", 
              "Ostrya", "Acer", "Acer", "Quercus", "Rubus", "Osmunda", "Fagus", "Prunus", 
              "Helianthus", "Solidago", "Agonis"),
    family = c("Myricaceae", "Sapindaceae", "Betulaceae", "Betulaceae", "Asteraceae", 
               "Salicaceae", "Salicaceae", "Rosaceae", "Osmundaceae", "Poaceae", "Poaceae", 
               "Sapindaceae", "Rosaceae", "Rosaceae", "Betulaceae", "Sapindaceae", 
               "Sapindaceae", "Fagaceae", "Rosaceae", "Osmundaceae", "Fagaceae", 
               "Rosaceae", "Asteraceae", "Asteraceae", "Myrtaceae")
  )
  
  # Prepare phylogenetic tree
  phylo <- phylo.maker(species_data, scenarios = "S3")
  phylo2 <- phylo.maker(species_data, scenarios = "S3")
  
  # Check for presence of species in phylo.maker
  phylo$scenario.3$tip.label
  
  # Resolve multichotomies and check tree properties
  #phylo <- multi2di(hypotheses$scenario.3)
  #phylo$scenario.3edge.length[phylo$edge.length <= 0] <- 1e-9
  
  # Compute the pairwise matrix of phylogenetic distances
  distance_matrix <- cophenetic.phylo(phylo2$scenario.3)
  
  # Return the distance matrix
  #return(as.matrix(distance_matrix))
#}
  
  write.csv(distance_matrix, "PD-s25-scenario3.csv")

  
### Print
  
  # Load required libraries
  library(ggplot2)
  library(reshape2)
  library(viridisLite)
  
  # Convert the distance matrix to a long-format data frame
  dist_df_long <- melt(as.matrix(distance_matrix))
  colnames(dist_df_long) <- c("species1", "species2", "distance")
  
  # Define the distance threshold to set break in legend
  threshold <- 300
  
  # Create a custom color palette
  custom_palette <- c(turbo(100), "gray90")  # 100 colors from turbo, plus white
  
  # Create a new fill variable in dist_df_long
  dist_df_long <- dist_df_long %>%
    mutate(
      fill_value = ifelse(distance > threshold, NA, distance),  # Set outliers to NA
      fill_color = ifelse(is.na(fill_value), "gray90", turbo(100)[as.numeric(cut(fill_value, breaks = 100))])
    )
  
#### PLOT
  ggplot(dist_df_long, aes(x = species1, y = species2, fill = fill_color)) +  # Use fill_color for the fill aesthetic
    geom_tile(color = "white") +  # Keep the tile borders white
    scale_fill_identity(  # Use scale_fill_identity to use the fill_color directly
      na.value = "gray90",  # Assign gray90 for NA values (outliers)
      name = "Phylogenetic\nDistance"  # Modify the title if needed
    ) + 
    coord_fixed() + 
    labs( 
      x = "Species", 
      y = "Species", 
      title = "Phylogenetic Distance Matrix" 
    ) + 
    theme_minimal(base_size = 10) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    annotate("text", x = Inf, y = Inf,  
             label = "Osmunda to rest: 781 M years divergent",  
             hjust = 1.15, vjust = 1.1, size = 3.5)
  
  
  
##### Antonio scripts
phylogenetic_diversity <- function(fia, 
                                   SPCD, 
                                   outpath, 
                                   progress,
                                   threads) {
  
  ###---------------------------------------------------------------------------
  # Merge inventory with taxonomic information
  fia_summary <- merge(fia, SPCD[, c("SPCD", "Taxon2", "Genus", "Family")], 
                       by = "SPCD", all.x = TRUE, all.y = FALSE)
  
  ### Prepare phylo  -----------------------------------------------------------
  # Summary of species and families
  species <- unique(fia_summary[, c("Taxon2", "Genus", "Family")])
  colnames(species) <- c("species", "genus", "family")
  species <- na.exclude(species)
  hypotheses <- phylo.maker(species, scenarios = "S3")
  
  # Resolve Multichotomies
  phylo <- multi2di(hypotheses$scenario.3)
  is.binary.phylo(phylo) #Test for Binary Tree
  is.ultrametric(phylo) #Test if a Tree is Ultrametric
  
  # Transform phylo
  tol = 1e-9
  phylo$edge.length[phylo$edge.length <= 0] <- tol
  is.ultrametric(phylo)
  #phylo$edge.length <- (log10(phylo$edge.length)+(min(log10(phylo$edge.length))*-1))
  
  #plot(phylo, type = "fan")

  ### Community matrix ---------------------------------------------------------
  
  # dcast
  master_matrix <- sample2matrix(fia_summary[, c("PLT_CN", "VOL_TPA", "Taxon2")])
  
  # Standardization of communities
  #master_matrix <- decostand(master_matrix, method = "total")
  #master_matrix <- round(master_matrix, 4)
  
  ### Get phylogenetic beta diversity ------------------------------------------
  
  #Get match
  matched <- picante::match.phylo.comm(phy = phylo, comm = master_matrix)
  
  cat(paste0("Number of communities:", nrow(master_matrix)))
  
  # Export
  saveRDS(matched$comm, paste0(out_path, "/master_matched_comm_phylo_VOL.rds"))
  saveRDS(matched$phy, paste0(out_path, "/master_matched_phylo.rds"))
  
  # Compute functional beta diversity
  beta <- beta_rcpp(comm = as.matrix(matched$comm), 
                    tree = matched$phy, 
                    taxonomic = FALSE,
                    threads = threads,
                    progress = progress)
  
  saveRDS(beta, paste0(out_path, "/beta_phylo_VOL.rds"))
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "E:/Projects/MMDB"
#root_path <- "/media/antonio/Work/MMDB"
root_path <- "/panfs/jay/groups/17/cavender/guzman/MMDB"

fia <- readRDS(paste0(root_path, "/03_reshape/fia_reshape.rds"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD_meta.csv"), na.strings="NA")
out_path <- paste0(root_path, "/03_diversity")
progress <- TRUE
threads <- 120

phylogenetic_diversity(fia, 
                       SPCD, 
                       outpath, 
                       progress,
                       threads)
