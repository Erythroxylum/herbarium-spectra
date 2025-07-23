#'------------------------------------------------------------------------------
#' @title-Get-Phylogenetic-Distances
#'------------------------------------------------------------------------------

#' @description Estimates distances from Smith & Brown tree using the 
#' V.PhyloMaker2 package
#' @return Return distance matrices as .rds files.


#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(phytools)
library(V.PhyloMaker2)

#'------------------------------------------------------------------------------
#' Source code


#'------------------------------------------------------------------------------
#' Function

######## Time Tree 5 (https://timetree.org/)

#Deal with missing taxa: Take this tree and add Phragmites australis as sister to Phalaris arundinacea at 39.751078, and add Betula populifolia as sister to Betula_papyrifera at 39.647898 (phylomaker v3 results)
phylo <- bind.tip(phylo, "Phragmites_australis", where = which(phylo$tip.label == "Phalaris_arundinacea"), position = 39.751078)
phylo <- bind.tip(phylo, "Betula_populifolia", where = which(phylo$tip.label == "Betula_papyrifera"), position = 39.647898)
phylo <- ladderize(phylo, right = FALSE)
plot(phylo, show.tip.label = TRUE, cex = 0.8)
phylo$node.label <- seq_along(phylo$node.label)

write.tree(phylo, "herbarium-predictors-analysis/phylogram_pd_TimeTree5.tre")


######## phylomaker

  # Read the dataset from a CSV file
data <- read.csv("DMWhiteHUHspec1_sp25leaf560_ref5nm_450-2400.csv")

# Extract relevant columns and retain unique species
species_data <- unique(data[, c("scientificName", "Genus", "Family")])
#species_data$species <- gsub(" ", "_", species_data$scientificName)

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

## Write trees

# Add underscores to tip labels
phylo$scenario.3$tip.label <- gsub(" ", "_", phylo$scenario.3$tip.label)

# Save the phylogram with PD branch lengths
phylogram <- phylo$scenario.3

phylogram$tip.label <- gsub("Osmunda_claytoniana", "Claytosmunda_claytoniana", phylogram$tip.label)

# Save tree as RDS files
saveRDS(phylogram, file = "herbarium-predictors-analysis/phylogram_pd.rds")   # Phylogram with PD branch lengths



#######################################
### PLOTS

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
custom_palette <- c(turbo(100), "gray90")# 100 colors from turbo, plus white

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
  
