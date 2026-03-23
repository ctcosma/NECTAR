#### Read in data ####
# all predicted interactions
predicted_interactions <- read.csv("Data_Clean/Interactions/Predicted_Interactions.csv")

# Calscape plant names
calscape_plants_harmonized <-  read.csv("Data_Clean/Species_Traits_Attributes/Calscape_Plants_Harmonized.csv")

#### Filter plants for only those in calscape, assign calscape names back ####
# Get all unique harmonized genus_species names, with original calscape name to change back to
calscape_plants_unique <- calscape_plants_harmonized %>%
  dplyr::select(genus_species_harm, genus_species) %>%
  unique()

# Join
predicted_interactions <- left_join(predicted_interactions, calscape_plants_unique, by = join_by("lower" == "genus_species_harm"), relationship = "many-to-many")

# assign calscape name to lower
predicted_interactions$lower <- predicted_interactions$genus_species

# filter for columns we want
predicted_interactions <- dplyr::select(predicted_interactions, c(JEPCODE, taxon, interactionTypeName, higher, lower, categoricalInteraction))

# Remove duplicate rows
predicted_interactions <- unique(predicted_interactions)

# filter out NAs in lower
predicted_interactions <- predicted_interactions[!is.na(predicted_interactions$lower), ]
  
## Save
write.csv(predicted_interactions, "Data_Calscape/Interactions/calscape_interactions_allPredicted.csv")
