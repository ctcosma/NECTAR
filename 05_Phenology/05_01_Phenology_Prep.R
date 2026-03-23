#### Load Packages ####
library(dplyr)
library(bdc)
library(stringr)

#### Read in necessary data ####
# Checklists
checklists <- read.csv("Data_Clean/Species_Checklists/checklist_cleaned.csv")

# Read in gbif annotated records
gbif_annotated_plants <- read.csv("Data_Raw/Phenology/annotations_all2.csv")

#### Filter to bbox around CA for faster processing ####
gbif_annotated_plants <- gbif_annotated_plants %>%
  filter(latitude > 32.5,
         latitude < 42.1,
         longitude > -124.5,
         longitude < 114.1)

#### Harmonize gbif_plant_annotations names ####
source("Code/99_Supporting/name_harmonization.R")

# 1. Parse names
phenobase_names <- unique(gbif_annotated_plants$scientific_name)
names_phenobase_parsed <- bdc_clean_names(sci_names = phenobase_names, save_outputs = FALSE)

# 2. Harmonize names
names_harm_phenobase <- harmonize_names(
  names = unique(names_phenobase_parsed$names_clean),
  higher_tax = "Plantae",
  names_file = "Temp/names_phenobase.tsv",
  names_harm_file = "Temp/names_harmonized_phenobase.tsv"
)

# 3. Apply Jepson taxonomy
# Note: some plants not appearing in CA will have names_harm == NA. This is okay,
# since our occurrence dataset only includes names in Jepson anyway.
names_harm_phenobase <- apply_jepson_taxonomy(names_harm_phenobase)

# Merge names with species lists
phenobase_harmonized <- gbif_annotated_plants %>%
  left_join(names_phenobase_parsed, by = join_by("scientific_name" == "scientificName")) %>%
  left_join(names_harm_phenobase, by = c("names_clean" = "ScientificName")) %>%
  unique()

#### Cleaning phenobase records ####
# Filter for only flowering records
phenobase_harmonized_clean <- phenobase_harmonized %>%
  filter(trait == "flower")

# Put in same format as occurrence data
phenobase_harmonized_clean <- phenobase_harmonized_clean %>%
  dplyr::select(latitude, longitude, names_harm, year, day_of_year, DataSourceTitle, ClassificationPath) %>%
  rename(scientificName = names_harm,
         dataSourceTitle = DataSourceTitle,
         classificationPath = ClassificationPath,
         yearClean = year,
         DOY = day_of_year,
         decimalLatitude = latitude,
         decimalLongitude = longitude) %>%
  mutate(taxon = "plants",
         phenology = "Flowering",
         phenologyCode = as.character("F"),
         databaseID = "gbif_annotations")

phenobase_harmonized_superclean <- phenobase_harmonized_clean %>%
  mutate(genus_species = str_extract(scientificName, "^[A-Z][a-z]*(?:\\s[a-z]+)?")) %>%
  filter(yearClean <= 2025,
         yearClean >= 1900,
         !is.na(decimalLatitude),
         !is.na(decimalLongitude),
         genus_species %in% checklists$genus_species)

# De-dupe
phenobase_harmonized_deduped <- phenobase_harmonized_superclean %>%
  group_by(genus_species, decimalLatitude, decimalLongitude, yearClean, DOY) %>%
  slice_head(n = 1) %>%
  ungroup()

# Write as csv
write.csv(phenobase_harmonized_deduped, "Data_Clean/Occurrences/gbif_plant_annotations_clean.csv")
