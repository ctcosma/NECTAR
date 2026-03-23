#### Load Packages ####
library(dplyr)

#### For Calscape Pheno Display ####
## read in data
# Phenometrics
phen_best <- read.csv("Data_Clean/Phenology/phenology_best_available.csv")
phen_simple <- read.csv("Data_Clean/Phenology/phenology_simple.csv")
phen_type2 <- read.csv("Data_Clean/Phenology/phenology_type2.csv")
phen_type3 <- read.csv("Data_Clean/Phenology/phenology_type3.csv")

# checklists
checklists <- read.csv("Data_Clean/Species_Checklists/checklist_cleaned.csv")

# # read in predicted interactions, use to remove pollinators not modelled
# predicted_interactions <- read.csv("Data_Clean/Interactions/Predicted_Interactions.csv")

## supporting vectors
month_names <- c("January", "Febuary", "March", "April", "May", "June", 
                 "July", "August", "September", "October", "November", "December")
month_breaks <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)

## Reformat phenology
# remove NA
phen_best <- na.omit(phen_best[, c("genus_species", "JEPCODE", "phen_start", "phen_end")])
phen_simple <- na.omit(phen_simple[, c("genus_species", "phen_start", "phen_end")])
phen_type2 <- na.omit(phen_type2[, c("genus_species", "phen_start", "phen_end")])
phen_type3 <- na.omit(phen_type3[, c("genus_species", "phen_start", "phen_end")])

# Merge simple, type2, and type3 (all are at same stat-level resolution)
phen_simple <- rbind(phen_simple, phen_type2, phen_type3)

# constrain to [0, 364]
phen_simple$phen_start[phen_simple$phen_start > 364] <- 364
phen_simple$phen_end[phen_simple$phen_end > 364] <- 364
phen_best$phen_start[phen_best$phen_start > 364] <- 364
phen_best$phen_end[phen_best$phen_end > 364] <- 364

# add taxa column
phen_simple <- inner_join(phen_simple,
                          unique(select(checklists, genus_species, taxon)), 
                          by = "genus_species")
phen_best <- inner_join(phen_best, 
                        unique(select(checklists, genus_species, taxon)), 
                        by = "genus_species")

## Pheno reformatting
get_month <- function(days) {
  findInterval(days, month_breaks)
}

phenomatrix <- function(phenoframe_onset, phenoframe_offset, month_names){
  # Create onset matrix (1s from start month to December)
  onset_matrix <- matrix(0, nrow = length(phenoframe_onset), ncol = 12)
  start_months <- get_month(phenoframe_onset)
  for(i in 1:nrow(onset_matrix)) {
    onset_matrix[i, start_months[i]:12] <- 1
  }
  
  # Create end matrix (1s from January to end month)
  end_matrix <- matrix(0, nrow = length(phenoframe_offset), ncol = 12)
  end_months <- get_month(phenoframe_offset)
  for(i in 1:nrow(end_matrix)) {
    end_matrix[i, 1:end_months[i]] <- 1
  }
  
  # Multiply matrices element-wise to get overlap
  phenomatrix <- onset_matrix * end_matrix
  colnames(phenomatrix) <- month_names
  
  return(phenomatrix)
}

# Calscape phenoframes
phen_simple_mat <- phenomatrix(phen_simple$phen_start, phen_simple$phen_end, month_names)
phen_simple_calscape <- cbind(phen_simple, phen_simple_mat)

phen_best_mat <- phenomatrix(phen_best$phen_start, phen_best$phen_end, month_names)
phen_best_calscape <- cbind(phen_best, phen_best_mat)

# Save
write.csv(phen_simple_calscape, "Data_Calscape/Phenology/calscape_phenology_stateScale.csv")
write.csv(phen_best_calscape, "Data_Calscape/Phenology/calscape_phenology_jepcodeScale.csv")
