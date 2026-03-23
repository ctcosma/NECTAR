#### Load packages ####
library(dplyr)
library(tibble)
library(data.table)
library(sf)
library(tigris)
library(CoordinateCleaner)
library(lubridate)
library(bdc)
library(BeeBDC)
library(stringr)

#### Read in harmonized occurrence database file ####
database <- fread(
  "Data_Raw/Formatted/database.csv",
  encoding = "UTF-8",      # ensure consistent encoding
  showProgress = TRUE,     # progress bar
  nThread = parallel::detectCores() - 1  # use all but one CPU core
)

#### Step 0: Manual Column Cleaning/Manipulation ####
# for all cch2 data downloaded using phenology flag, mark as flowering, if not already
database <- database %>%
  mutate(
    Phenology = if_else(grepl("cch2_pheno", database_id, ignore.case = TRUE), "Flowering", Phenology),
    Phenology.Code = if_else(grepl("cch2_pheno", database_id, ignore.case = TRUE), "F", Phenology.Code)
  )

#### Step 1: Basic column + spatial data ####
## column checks
# Records missing species names
check_pf <-
  bdc_scientificName_empty(
    data = database,
    sci_name = "scientificName")

# Records lacking information on geographic coordinates
check_pf <- bdc_coordinates_empty(
  data = check_pf,
  lat = "decimalLatitude",
  lon = "decimalLongitude")

# Records with out-of-range coordinates
check_pf <- bdc_coordinates_outOfRange(
  data = check_pf,
  lat = "decimalLatitude",
  lon = "decimalLongitude")

# Records from poor sources
check_pf <- bdc_basisOfRecords_notStandard(
  data = check_pf,
  basisOfRecord = "basisOfRecord",
  names_to_keep = c("Event","HUMAN_OBSERVATION", "HumanObservation", 
                    "LIVING_SPECIMEN", "LivingSpecimen", "MACHINE_OBSERVATION", 
                    "MachineObservation", "MATERIAL_SAMPLE", "None", "O", 
                    "Occurrence", "MaterialSample", "OBSERVATION", 
                    "Pinned Specimen", "Photograph", "Preserved Specimen", 
                    "PRESERVED_SPECIMEN", "preservedspecimen Specimen", 
                    "Preservedspecimen", "PreservedSpecimen", 
                    "preservedspecimen", "S", "Specimen", "Taxon", 
                    "UNKNOWN", "", NA))

# Since all records in US, to see if any lat/lon transposed
# Ensure no negative numbers in latitude and no positives in longitude
check_pf$.coordinates_transposed = ifelse(check_pf$decimalLatitude<0, FALSE, TRUE)

## Spatial cleaning

# Get map sf object for CA + 100km buffer
options(tigris_use_cache = TRUE)

# 1) Get California boundary (as sf polygon)
ca <- states(cb = FALSE) %>%
  filter(STUSPS == "CA") %>%
  st_make_valid()

# 2) Reproject to California Albers (EPSG:3310) for accurate distance buffering
ca_albers <- st_transform(ca, 3310)

# 3) Create 100 km buffer
ca_buffer_100km <- st_buffer(ca_albers, dist = 100000)

# 4) Transform back to WGS84 (EPSG:4326) for lat/long use
CA_boundary <- st_transform(ca_buffer_100km, 4326)

#fix invalid geometries
CA_boundary <- st_make_valid(CA_boundary)

# Check if points are in California + 100km buffer (region of interest)
check_pf <- st_as_sf(check_pf, 
                     coords = c("decimalLongitude", "decimalLatitude"), 
                     crs = st_crs(CA_boundary), 
                     remove = FALSE, 
                     na.fail = FALSE) %>%
  # Use st_intersects() to check for overlap between points and polygon
  # If no intersection, set .overlaps_California to FALSE
  mutate(.overlaps_California = lengths(st_intersects(., CA_boundary)) > 0) %>%
  # Remove geometry to turn back into standard dataframe
  st_drop_geometry()

# Create summary column
check_pf <- bdc_summary_col(data = check_pf)

# Create report
report <-
  bdc_create_report(data = check_pf,
                    database_id = "database_id",
                    workflow_step = "prefilter",
                    save_report = FALSE)
report

## Filter out flagged records
# Save flagged records in a separate dataset
flagged_records <- check_pf %>%
  dplyr::filter(.summary == FALSE) %>%
  dplyr::select(-.summary)

# Filter out all flagged records
merged_occurrences <-
  check_pf %>%
  dplyr::filter(.summary == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")

#### Part 2: Cleaning and Harmonizing Scientific Names ####
source("Code/99_Supporting/name_harmonization.R")

# Be sure to first install gnparser from https://github.com/gnames/gnparser
parse_names <-
  bdc_clean_names(sci_names = merged_occurrences$scientificName, save_outputs = FALSE)

# Examine errors
wrong_names <- filter(parse_names, is.na(names_clean) | quality == 0)

# Replace cleaned names in merged_occurrences file
parse_names <-
  parse_names %>%
  dplyr::select(.uncer_terms, names_clean)

merged_occurrences <- dplyr::bind_cols(merged_occurrences, parse_names)

### Name harmonization
# 1. Separate animals and plants due to possible duplication of names across Animalia and Plantae
merged_occurrences_plants <- merged_occurrences %>%
  filter(grepl("calflora|cch2", database_id, ignore.case = TRUE))
merged_occurrences_animals <- merged_occurrences %>%
  filter(!grepl("calflora|cch2", database_id, ignore.case = TRUE))

names_plants <- merged_occurrences_plants %>% pull(names_clean) %>% unique()
names_animals <- merged_occurrences_animals %>% pull(names_clean) %>% unique()

# 2. Harmonize names
names_harm_plants <- harmonize_names(
  names = names_plants,
  higher_tax = "Plantae",
  names_file = "Temp/names_occurrences_plants.tsv",
  names_harm_file = "Temp/names_harmonized_occurrences_plants.tsv"
)

names_harm_animals <- harmonize_names(
  names = names_animals,
  higher_tax = "Animalia",
  names_file = "Temp/names_occurrences_animals.tsv",
  names_harm_file = "Temp/names_harmonized_occurrences_animals.tsv"
)

# 3. Apply Jepson taxonomy to plants
names_harm_plants <- apply_jepson_taxonomy(names_harm_plants)

# 4. Join back with occurrence data
merged_occurrences_harmonized_plants <- left_join(merged_occurrences_plants,
                                                  names_harm_plants,
                                                  by = c("names_clean" = "ScientificName"))
merged_occurrences_harmonized_animals <- left_join(merged_occurrences_animals,
                                                   names_harm_animals,
                                                   by = c("names_clean" = "ScientificName"))

#### Part 3: Update Taxonomy columns ####
# add higher taxonomy with function
merged_occurrences_harmonized_plants <- update_taxonomy(merged_occurrences_harmonized_plants, tax = "plant")
merged_occurrences_harmonized_animals <- update_taxonomy(merged_occurrences_harmonized_animals, tax = "animal")

# combine
merged_occurrences_harmonized <- rbind(merged_occurrences_harmonized_plants, merged_occurrences_harmonized_animals)

# Add a taxa column that lists what specific taxa we are dealing with
merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  mutate(taxon = case_when(
    grepl("Plantae", kingdom, ignore.case = TRUE)  ~ "plants",
    grepl("trochilidae", family, ignore.case = TRUE)  ~ "hummingbirds",
    grepl("Apidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Andrenidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Colletidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Halictidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Megachilidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Mellittidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Syrphidae", family, ignore.case = TRUE)  ~ "hoverflies",
    grepl("Papilionidae", family, ignore.case = TRUE)  ~ "butterflies",
    grepl("Nymphalidae", family, ignore.case = TRUE)  ~ "butterflies",
    grepl("Pieridae", family, ignore.case = TRUE)  ~ "butterflies",
    grepl("Lycaenidae", family, ignore.case = TRUE)  ~ "butterflies",
    grepl("Riodinidae", family, ignore.case = TRUE)  ~ "butterflies",
    grepl("Hesperiidae", family, ignore.case = TRUE)  ~ "butterflies",
    grepl("Heliconiidae", family, ignore.case = TRUE)  ~ "butterflies",
    grepl("Lepidoptera", order, ignore.case = TRUE)  ~ "moths"
    ))

#examine errors
group_miss = filter(merged_occurrences_harmonized, is.na(taxon))

### Filter out flagged records and save checkpoint ###

# Add more flagged records to the flagged_records dataset
merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  # Create flag column for taxon group. FALSE if taxon group not in scope of project
  mutate(.correct_taxon = ifelse(!is.na(taxon), TRUE, FALSE)) %>%
  # Create flag column for harmonized name. remove if no resulting name
  mutate(.name_missing = ifelse(!is.na(names_harm), TRUE, FALSE))

# Filter out records with flag
incorrect_group <- merged_occurrences_harmonized %>%
  filter(.correct_taxon == FALSE,
         .name_missing == FALSE)

# Add the records with incorrect taxonomy to flagged records
flagged_records <- bind_rows(flagged_records, incorrect_group) %>%
  dplyr::select(-starts_with("."), starts_with(".")) # arrange the columns so flag columns are at the end!

# Remove any occurrences that don't fall into one of our taxon groups of interest
# These are going to be arthropods that are not bees, leps, or hoverflies
merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  filter(.correct_taxon == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")

# Save as csv as checkpoint
write.csv(merged_occurrences_harmonized, "Temp/merged_occurrences_harmonized_temp1.csv", row.names = F)
write.csv(flagged_records, "Temp/flagged_records_temp1.csv", row.names = F)

# Clean up working directory
rm(list = setdiff(ls(), c("merged_occurrences_harmonized", "flagged_records")))

#### Part 4: Spatial Cleaning ####
# # re-load files if needed
# merged_occurrences_harmonized = read.csv("Temp/merged_occurrences_harmonized_temp1.csv")
# flagged_records = read.csv("Temp/flagged_records_temp1.csv")

# Make coordinates as type numeric (instead of character)
merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  # make coordinates numeric instead of character
  mutate(decimalLatitude = as.numeric(decimalLatitude),
         decimalLongitude = as.numeric(decimalLongitude))

# Check other common spatial issues
check_space <- BeeBDC::jbd_coordinates_precision(data = merged_occurrences_harmonized, 
                                                 lon = "decimalLongitude", 
                                                 lat = "decimalLatitude", 
                                                 ndec = 1) %>%
  # Rename flag column to make more explicit
  rename(.coordinates_Precise = .rou)

# 1 decimal point is ~10km accuracy in CA, 2 is ~1km accuracy. Keeping to 1 for now, its not too far off from our SDM predictor grid resolution
# ^ This flagged nearly 2 million records when ndec=2

# # Removing this step - Most values are NA or suspect (VERY large, unlikely to be legitimate)
# # Filter or flag records with excessive coordinate uncertainty (keep NA values)
# check_space$.coordinates_uncertainty <-
#   is.na(check_space$coordinateUncertaintyInMeters) |
#   check_space$coordinateUncertaintyInMeters <= 10000
# 
# # Count records flagged as too uncertain (FALSE = failed the check)
# n_uncertain <- sum(!check_space$.coordinates_uncertainty)
# 
# # Calculate proportion
# prop_uncertain <- mean(!check_space$.coordinates_uncertainty)
# 
# # Display results
# cat("Records flagged as too uncertain:", n_uncertain, "\n")
# cat("Proportion flagged:", round(prop_uncertain * 100, 2), "%\n")

# Force countrycode to US so it does not flag anything because of this below
# (lots of records have no country code due to dataset it is coming from)
check_space$countryCode = "US"

#subset the dataset to run in chunks

#remove spent dataxset
merged_occurrences_harmonized = NULL

# Define chunk size
chunk_size <- 1e6  # 2 million rows

# Initialize final dataframe
check_space2 <- tibble()

# Total number of rows
total_rows <- nrow(check_space)

# Process in chunks
for (start_row in seq(1, total_rows, by = chunk_size)) {
  # Define the end row for the current chunk
  end_row <- min(start_row + chunk_size - 1, total_rows)
  
  # Subset the data for the current chunk
  chunk <- check_space[start_row:end_row, ]
  
  # Run the cleaning function on the current chunk
  chunk_cleaned <-
    CoordinateCleaner::clean_coordinates(
      x =  chunk,
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "names_harm",
      tests = c(
        "capitals",     # records within 0.5 km of capitals centroids
        "centroids",    # records within 1 km around country and province centroids
        "equal",      # records with equal coordinates
        "gbif",         # records within 1 km of GBIF headquarters. (says 1 degree in package, but code says 1000 m)
        "institutions", # records within 100m of zoo and herbaria
        "zeros"       # records with coordinates 0,0
        # "seas"        # Not flagged as this should be flagged by coordinate country inconsistent
      ),
      capitals_rad = 1000,
      centroids_rad = 500,
      centroids_detail = "both", # test both country and province centroids
      inst_rad = 100, # remove zoo and herbaria within 100m
      range_rad = 0,
      zeros_rad = 0.5,
      capitals_ref = NULL,
      centroids_ref = NULL,
      country_ref = NULL,
      country_refcol = "countryCode",
      inst_ref = NULL,
      range_ref = NULL,
      # seas_scale = 50,
      value = "spatialvalid" # result of tests are appended in separate columns
    ) %>%
    # Remove duplicate .summary column that can be replaced later and turn into a tibble
    dplyr::select(!tidyselect::starts_with(".summary")) %>%
    dplyr::tibble()
  
  # Append the cleaned chunk to the final dataframe
  check_space2 <- bind_rows(check_space2, chunk_cleaned)
  
  # Optional: print progress
  print(paste("Processed rows", start_row, "to", end_row))
}

### Filter out flagged records and save checkpoint ###
# Summarize flags
check_space_flagSummary <- BeeBDC::summaryFun(data = check_space2, 
                                              dontFilterThese = c(".uncer_terms"), 
                                              removeFilterColumns = FALSE,
                                              filterClean = FALSE)
  
# Convert data type of coordinates back to character for consistency with other objects
check_space_flagSummary <- check_space_flagSummary %>%
  mutate(decimalLatitude = as.character(decimalLatitude),
         decimalLongitude = as.character(decimalLongitude),
         individualCount = as.character(individualCount),
         endDayOfYear = as.character(endDayOfYear))

flagged_records <- flagged_records %>%
  mutate(decimalLatitude = as.character(decimalLatitude),
         decimalLongitude = as.character(decimalLongitude),
         individualCount = as.character(individualCount),
         endDayOfYear = as.character(endDayOfYear))

# Filter out records with flag
spatial_issues <- check_space_flagSummary %>%
  filter(.summary == FALSE)

# Add the spatial flagged records to flagged records
flagged_records <- bind_rows(flagged_records, spatial_issues) %>%
  dplyr::select(-starts_with("."), starts_with(".")) %>% # arrange the columns so flag columns are at the end!
  dplyr::select(-.summary)

# Remove the flagged rows
space_clean <- check_space_flagSummary %>%
  filter(.summary == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")

# Clean up working directory
rm(list = setdiff(ls(), c("space_clean", "flagged_records")))

# Save a copy as check point
write.csv(space_clean, "Temp/space_clean.csv", row.names = F)
write.csv(flagged_records, "Temp/flagged_records_temp2.csv", row.names = F)

#### Part 5: Temporal Cleaning #### 
# Reload data if needed 
# space_clean = read.csv("Temp/space_clean.csv")
# flagged_records = read.csv("Temp/flagged_records_temp2.csv")

#remove the mpg data, not considering those dates
mpg_rows <- space_clean[grepl("mpg", space_clean$database_id), ]
space_clean <- space_clean[!grepl("mpg", space_clean$database_id), ]

# Date cleaning - get date information from multiple columns and clean
occurrences_time <- space_clean %>%
  mutate( # First, add dates to clean columns using the date column, which has variable format depending on data source
    dateClean = case_when(
      grepl("\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z", eventDate) ~ as.Date(eventDate, format = "%Y-%m-%dT%H:%M:%OSZ"),
      grepl("\\d{4}-\\d{2}-\\d{2}", eventDate) ~ as.Date(eventDate, format = "%Y-%m-%d"),
      TRUE ~ as.Date(NA)
    ),
    yearClean = as.numeric(lubridate::year(dateClean)),
    monthClean = as.numeric(lubridate::month(dateClean)),
    dayClean = as.numeric(lubridate::day(dateClean))
  ) %>%
  mutate( # Then replace any NAs with values from the year/month/day column if values exist there
    yearClean = ifelse(is.na(yearClean) & !is.na(year), as.numeric(year), yearClean),
    monthClean = ifelse(is.na(monthClean) & !is.na(month), as.numeric(month), monthClean),
    dayClean = ifelse(is.na(dayClean) & !is.na(day), as.numeric(day), dayClean)
  ) %>%
  dplyr::select(-dateClean)

# Check that the warning message is not problematic

### Dataset Cleaning ###
# Filter out based on cleaner dates
occurrences_time_cleaned <- occurrences_time %>%
  mutate(.valid_date = case_when(
    yearClean >= 1900 & yearClean <= 2025 ~ TRUE,
    TRUE ~ FALSE
  ))

# Filter out records with flag
temporal_issues <- occurrences_time_cleaned %>%
  filter(.valid_date == FALSE)

#examine flagged_records to see if time cleaning is removing valid rows
time = temporal_issues

calflora <- time[grepl("calflora", time$database_id), ]
calflora = calflora %>%
  dplyr::select("eventDate", "day", "month", "year", "dayClean", "monthClean", "yearClean")

cch2 <- time[grepl("cch2", time$database_id), ]
cch2 = cch2 %>%
  dplyr::select("eventDate", "day", "month", "year", "dayClean", "monthClean", "yearClean")

gbif <- time[grepl("gbif", time$database_id), ]
gbif = gbif %>%
  dplyr::select("eventDate", "day", "month", "year", "dayClean", "monthClean", "yearClean")

ecdysis <- time[grepl("ecdysis", time$database_id), ]
ecdysis = ecdysis %>%
  dplyr::select("eventDate", "day", "month", "year", "dayClean", "monthClean", "yearClean")

chesshire <- time[grepl("chesshire", time$database_id), ]
chesshire = chesshire %>%
  dplyr::select("eventDate", "day", "month", "year", "dayClean", "monthClean", "yearClean")

idigbio <- time[grepl("idigbio", time$database_id), ]
idigbio = idigbio %>%
  dplyr::select("eventDate", "day", "month", "year", "dayClean", "monthClean", "yearClean")

#convert column to character
flagged_records$individualCount = as.character(flagged_records$individualCount)

# Add the records with incorrect taxonomy to flagged records
flagged_records <- bind_rows(flagged_records, temporal_issues) %>%
  dplyr::select(-starts_with("."), starts_with(".")) # arrange the columns so flag columns are at the end!

# Remove the flagged rows
time_clean <- occurrences_time_cleaned %>%
  filter(.valid_date == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")

# combine again with mpg_rows
time_clean = bind_rows(time_clean, mpg_rows)

# Clean up working directory
rm(list = setdiff(ls(), c("time_clean", "flagged_records")))

###Save checkpoint###
write.csv(time_clean, "Temp/time_clean.csv", row.names = F)
write.csv(flagged_records, "Temp/flagged_records_temp3.csv", row.names = F)

#### Part 6: De-duplication ####
# # Reload data if needed 
time_clean = read.csv("Temp/time_clean.csv")
flagged_records = read.csv("Temp/flagged_records_temp3.csv")

# Fix database_id column
time_clean <- time_clean %>%
  mutate(database_id = sub("_[^_]*$", "", database_id)) %>%
  mutate(database_id = sub("_formatted", "", database_id))

#fix problematic columns
time_clean$individualCount <- as.character(time_clean$individualCount)
time_clean$endDayOfYear = as.character(time_clean$endDayOfYear)

## Good resource on identifying duplicates:
# https://discourse.gbif.org/t/duplicate-occurrence-records/3735

## 3 general types
# exact duplicates have the same entries in all fields, including unique identifiers like occurrenceID.
# strict duplicates have the same entries in all fields except those which give the record a unique ID, like the occurrenceID field.
# relaxed duplicates have the same entries in key occurrence fields, such as taxon, coordinates, date and collector/observer.

#Chris's updated duplicate code as of Nov 22, 2024

#remove the mpg data, we are assuming those are already de-duplicated
mpg_rows <- time_clean[grepl("mpg", time_clean$database_id), ]

time_clean <- time_clean[!grepl("mpg", time_clean$database_id), ]

#Fix the phenology code issue - data type where F (for "Flowering") read as False (Boolean)
time_clean <- time_clean %>%
  mutate(Phenology.Code = if_else(Phenology.Code == FALSE, "F", Phenology.Code))

# Split the dataset into each taxonomic group
#Compile a list of dataframes for de-duplication

df_list <- list(
  
  #plant data for phenometrics
  plants_pheno = filter(time_clean, taxon == "plants" & (!(is.na(Phenology.Code)) | Phenology.Code != "")),
  
  #Calflora plants for SDMs
  # plants = filter(time_clean, database_id == "calflora"),
  
  # All plants for SDMs
  all_plants = filter(time_clean, database_id %in% c("calflora", "cch2")),
  
  # Chesshire for bees
  bees_chesshire = filter(time_clean, database_id == "chesshire"),
  
  # # All bees
  bees_all = filter(time_clean, taxon == "bees"),
  
  # Hoverflies
  hoverflies = filter(time_clean, taxon == "hoverflies"),
  
  # Hummingbirds
  hummingbirds = filter(time_clean, taxon == "hummingbirds"),
  
  # Leps
  leps = filter(time_clean, taxon %in% c("butterflies", "moths"))
  
)

# Remove all unneeded datasets
rm(list = setdiff(ls(), c("df_list", "flagged_records", "mpg_rows")))

##### Identify and filter out duplicates #####

# initiate blank df to fill
clean_occurrences = data.frame()

#fix column
flagged_records$individualCount = as.character(flagged_records$individualCount)
flagged_records$endDayOfYear = as.character(flagged_records$endDayOfYear)

#loop through each for de-duplication within each

for (i in 1:length(df_list)) {
  
  #identify df
  df = df_list[[i]]
  
  # Convert the dataset to a data.table object because of memory/computing limits
  df <- as.data.table(df)
  
  # Determine the duplicate types
  
  #Type 1: .dupExact = rows that are duplicates in everything (except "database_id", which is a column we added)
  #Type 2: .dupStrict = rows that are duplicates in everything except "database_id","occurrenceID"
  #Type 3: .dupRelaxed = rows that are duplicates based only on "names_harm", "decimalLatitude", "decimalLongitude", "yearClean", "monthClean", "dayClean", "Phenology", "Phenology.Code"
  df[, `:=`(
    .dupExact = duplicated(.SD, by = names(df)[!names(df) %in% "database_id"]),
    
    .dupStrict = duplicated(.SD, by = names(df)[!names(df) %in% c("database_id", "occurrenceID")]),
    
    dupRelaxed = duplicated(.SD, by = c("names_harm", "decimalLatitude", "decimalLongitude", 
                                        "yearClean", "monthClean", "dayClean", 
                                        "Phenology", "Phenology.Code"))
  )]
  
  # Filter out records with flag (keeping first occurrence of duplicates)
  # Not removing dupRelaxed for now
  duped_remove <- df %>%
    filter(.dupExact == TRUE | .dupStrict == TRUE)
  
  # Add the records with duplicate flag to flagged records
  flagged_records <- bind_rows(flagged_records, duped_remove) %>%
    dplyr::select(-starts_with("."), starts_with(".")) # arrange the columns so flag columns are at the end!
  
  # Remove the flagged rows (keeping first occurrence)
  df <- df %>%
    filter(.dupExact == FALSE & .dupStrict == FALSE) %>%
    bdc_filter_out_flags(data = ., col_to_remove = "all")
  
  # Add in a category name for repeat datasets (we are doing separate versions for SDMs vs checklists)
  df$category = names(df_list)[i]
  
  #Combine to single dataframe
  clean_occurrences = bind_rows(clean_occurrences, df)
  
}

#remove spent dataset
df_list = NULL

#combine back with mpg_rows
clean_occurrences = bind_rows(clean_occurrences, mpg_rows)

#### Part 7: Removing absence records ####

# identify absences to flagged_records?
clean_occurrences <- clean_occurrences %>%
  mutate(.presence_record = case_when(
    individualCount != 0 | is.na(individualCount) ~ TRUE,
    TRUE ~ FALSE
  ))

# filter absences
absences <- clean_occurrences %>%
  filter(.presence_record == F)

# Add the absence records to flagged
flagged_records <- bind_rows(flagged_records, absences) %>%
  dplyr::select(-starts_with("."), starts_with(".")) # arrange the columns so flag columns are at the end!

# Remove the flagged rows
clean_occurrences <- clean_occurrences %>%
  filter(.presence_record == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")


#### Part 8: Column filtering ####
#Fix mpg rows
clean_occurrences <- clean_occurrences %>%
  mutate(category = if_else(database_id == "mpg", "leps", category))

# filter out if only genus
clean_occurrences <- clean_occurrences %>%
  # how many parts of scientific name?
  mutate(num_parts_name = sapply(strsplit(trimws(names_harm), "\\s+"), length)) %>%
  # filter out if only 1 (genus/family)
  filter(num_parts_name > 1) %>%
  # remove previously created column
  dplyr::select(-num_parts_name)

# Create genus_species column
clean_occurrences <- clean_occurrences %>%
  mutate(genus_species = sapply(strsplit(names_harm, "\\s+"), function(x) paste(x[1:2], collapse = " ")))

# Only select what columns we need
occurrences_clean_niceColumns <- clean_occurrences %>%
  dplyr::select(category, database_id, taxon, names_harm, genus_species, decimalLatitude, decimalLongitude, yearClean, monthClean, dayClean, sex, lifeStage, Phenology, Phenology.Code, kingdom, phylum, class, order, superfamily, family, tribe, subtribe, genus, specificEpithet, infraspecificEpithet) %>%
  rename(scientificName = names_harm,
         databaseID = database_id,
         phenology = Phenology,
         phenologyCode = Phenology.Code)

##### Save Outputs #####
# Save occurrences
write.csv(occurrences_clean_niceColumns, "Temp/occurrences_clean.csv", row.names = F)

# Save occurrences - all columns
write.csv(clean_occurrences, "Temp/occurrences_clean_allColumns.csv", row.names = F)

# Save flagged records
write.csv(flagged_records, "Temp/occurrences_flagged.csv", row.names = F)
