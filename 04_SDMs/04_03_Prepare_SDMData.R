##### Prepare Data for SDMs #####
## Load data

library(data.table)
library(dplyr)
library(sf)

# setwd()

# Load full cleaned occurrence data
occurrences_clean = read.csv("Data_Clean/Occurrences/occurrences_clean.csv")

#prepare dataset
sdm_data = occurrences_clean %>%
  filter(!(category == "plants_pheno")) %>%  #remove pheno plants
  dplyr::select(taxon, genus_species, decimalLongitude, decimalLatitude, scientificName) %>% #select needed columns
  rename(
    taxon = taxon,
    species = genus_species,
    longitude = decimalLongitude,
    latitude = decimalLatitude,
    species_full = scientificName
  ) #rename columns

## Load California boundary (e.g., from rnaturalearth)
# You can substitute your own shapefile if preferred:
# ca <- st_read("Data_Spatial/CA_Boundary/CA_state_boundary.shp")
ca <- rnaturalearth::ne_states(country = "united states of america", returnclass = "sf") %>%
  filter(name == "California") %>%
  st_transform(4326)

## Convert occurrences to sf points
sdm_sf <- st_as_sf(sdm_data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

## Spatial filter: keep points within California polygon (takes a very long time with large datasets)
sdm_sf_in_CA <- sdm_sf[ca, , op = st_within]

## Identify species with ≥1 point inside CA
species_with_CA <- unique(sdm_sf_in_CA$species)

## Filter the full dataset to keep only those species
sdm_data_final <- sdm_data %>%
  filter(species %in% species_with_CA)

##### Save Outputs #####
write.csv(sdm_data_final, "Data_Clean/SDMs/SDMs_Inputs/sdm_data.csv", row.names = F)

## split into chunks by species and save
unique_species <- unique(sdm_data_final$species)  

# Assign each species to a group (1-11)
species_groups <- data.frame(
  species = unique_species,
  group = rep(1:11, length.out = length(unique_species))
)

# Join and split
sdm_data_with_groups <- sdm_data_final %>%
  left_join(species_groups, by = "species") 

# Split into list of 11 dataframes
sdm_data_list <- split(sdm_data_with_groups %>% dplyr::select(-group), 
                       sdm_data_with_groups$group)

# Save each df
for(i in 1:11) {
  write.csv(sdm_data_list[[i]], 
            file = paste0("Data_Clean/SDMs/SDMs_Inputs/occ_split/sdm_data_", i, ".csv"), 
            row.names = FALSE)
}
