#### Load packages ####
library(dplyr)
library(stringr)
library(sf)
library(terra)
library(data.table)

#### Read in data ####
# California ecoregions
ecoregions <- read_sf("Data_Raw/Spatial/Jepson_Ecoregions/jepson.shp") %>%
  st_make_valid()

## Checklists
# all checklists combined
checklist_cleaned <- fread("Data_Clean/Species_Checklists/checklist_cleaned.csv")

# sp x region checklist
sp_by_region <- read.csv("Data_Clean/Species_Checklists/ecoregion_checklist.csv")

# Plant checklist
plant_checklist <- checklist_cleaned[checklist_cleaned$taxon == "plants", ] %>%
  dplyr::select(genus_species) %>%
  mutate(genus = str_extract(genus_species, "[A-Z][a-z]*")) %>%
  unique()

# Pollinator checklists
pollinator_checklists <- checklist_cleaned[checklist_cleaned$taxon != "plants", ] %>%
  dplyr::select(genus_species) %>%
  mutate(genus = str_extract(genus_species, "[A-Z][a-z]*")) %>%
  unique()

# Occurrences
occurrences <- fread("Data_Clean/Occurrences/occurrences_clean.csv")

# Phenology
phenology <- fread("Data_Clean/Phenology/phenology_best_available.csv") %>% distinct()

## Interactions
interactions <- read.csv("Data_Clean/Interactions/ints_final_clean_slim.csv") %>%
  # remove hummingbird interactions from raw data
  filter(!(sourceTaxonGenusName %in% checklist_cleaned$genus[checklist_cleaned$taxon == "hummingbirds"]))
interactions_nectar <- filter(interactions, interactionTypeName == "visitsFlowersOf")
interactions_host   <- filter(interactions, interactionTypeName == "hasHost")

predicted_interactions <- fread("Data_Clean/Interactions/Predicted_Interactions.csv") %>%
  filter(sourceTaxonType != "hummingbirds")

#### Create spatial datasets ####
# Helper to convert occurrence/interaction df to spatial and join ecoregions
add_ecoregion <- function(df, lon_col = "decimalLongitude", lat_col = "decimalLatitude") {
  df %>%
    mutate(
      !!lon_col := as.numeric(.data[[lon_col]]),
      !!lat_col := as.numeric(.data[[lat_col]])
    ) %>%
    filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]])) %>%
    mutate(
      orig_longitude = .data[[lon_col]],
      orig_latitude  = .data[[lat_col]]
    ) %>%
    st_as_sf(coords = c(lon_col, lat_col), crs = 4326) %>%
    st_transform(crs(ecoregions)) %>%
    st_join(ecoregions, join = st_intersects) %>%
    mutate(
      !!lon_col := orig_longitude,
      !!lat_col := orig_latitude
    ) %>%
    dplyr::select(-orig_longitude, -orig_latitude)
}

occurrences_spatial  <- add_ecoregion(occurrences)
interactions_spatial <- add_ecoregion(interactions)

# Filter for CA native interactions
native_CA_interactions <- interactions_spatial %>%
  filter(sourceTaxonName_harm %in% pollinator_checklists$genus_species &
           targetTaxonName_harm %in% plant_checklist$genus_species &
           !is.na(JEPCODE))

predicted_interactions_onlySpatiallyExtrapolated <- predicted_interactions %>%
  filter(sourceTaxonName_harm %in% pollinator_checklists$genus_species &
           targetTaxonName_harm %in% plant_checklist$genus_species) %>%
  filter(paste(sourceTaxonName_harm,
               targetTaxonName_harm,
               interactionTypeName) %in%
           paste(interactions$sourceTaxonName_harm,
                 interactions$targetTaxonName_harm,
                 interactions$interactionTypeName))
