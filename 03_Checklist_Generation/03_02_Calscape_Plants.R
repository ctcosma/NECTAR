#### Load packages ####
library(dplyr)
library(bdc)

#### Load Data ####
# Calscape plant names
calscape_plants <- read.csv("Data_Raw/Species_Traits_Attributes/plants_table_5-7.csv")

#### Harmonize calscape plant names ####
# Load harmonization functions
source("Code/99_Supporting/name_harmonization.R")

# 1. Parse names
calscape_names <- unique(calscape_plants$species)
names_calscape_parsed <- bdc_clean_names(sci_names = calscape_names, save_outputs = FALSE)

# 2. Harmonize names
names_harm_calscape <- harmonize_names(
  names = unique(names_calscape_parsed$names_clean),
  higher_tax = "Plantae",
  names_file = "Temp/names_calscape.tsv",
  names_harm_file = "Temp/names_harmonized_calscape.tsv"
)

# 3. Apply Jepson taxonomy
# Note: some plants not appearing in CA will have names_harm == NA. This is okay,
# since our occurrence dataset only includes names in Jepson anyway.
names_harm_calscape <- apply_jepson_taxonomy(names_harm_calscape)

# Merge names with species lists
calscape_plants_harmonized <- calscape_plants %>%
  # first with bdc cleaned names
  left_join(names_calscape_parsed, by = join_by("species" == "scientificName")) %>%
  # then with harmonized names
  left_join(names_harm_calscape, by = c("names_clean" = "ScientificName")) %>%
  unique()

# Create genus_species column of harmonized name and original
calscape_plants_list <- calscape_plants_harmonized %>%
  mutate(genus_species = sapply(strsplit(species, "\\s+"), function(x) paste(x[1:2], collapse = " ")),
         genus_species_harm = sapply(strsplit(names_harm, "\\s+"), function(x) paste(x[1:2], collapse = " ")))

#### Edit this step if/when we want to treat ssp and vars differently
# Remove cultivars, default to species level
calscape_plants_list <- filter(calscape_plants_list, is_cultivar == 0)

# Keep only name columns
calscape_plants_list <- calscape_plants_list %>%
  dplyr::select(c("species", "genus_species", "names_harm", "genus_species_harm")) %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  unique()

#### Save plant supporting data with harmonized names ####
write.csv(calscape_plants_list, "Data_Clean/Species_Traits_Attributes/Calscape_Plants_Harmonized.csv", row.names = F)