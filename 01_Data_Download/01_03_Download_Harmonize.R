#### Load Packages ####
library(purrr)
library(readr)
library(sf)
library(tigris)
library(rmapshaper)
library(dplyr)
library(rgbif)
library(data.table)
library(ridigbio)

#### Criteria ####
## Sources for each taxa
# Plants: Calflora, CCH2
# Pollinators: GBIF, Ecydysis (formerly SCAN), Chesshire (bees), idigbio

## Download criteria: CA + 100 km buffer (except for CCH2 pheno)

#### Plants ####

##### CalFlora #####
### Manual Approach
# Go to https://www.calflora.org/entry/wsearch.html 
# Need a CalFlora Account to use the downloader
# Change "Native Status" to  "Native"
# Change to Wild status
# Change "Download Format" to "CSV"
# Change "Geometry" to "point"
# Specify columns to include: ID	Taxon	Latitude	Longitude	Date	Observer	Location Description	Habitat	Phenology	Phenology Code	Source	County	Accuracy: Square Meters	Error Radius (m)	Location Quality
# Click on each County from the list one at a time (download limit of 500,000 records is reached when doing "Any") then click "Search", then click "Download". Skip "- Bay Area -"
# Should be 58 total files

### Automated Approach

# 1) run script 01_01_Calflora.R (if not done already)

# 2) Read in all files, merge
calflora_files <- list.files("Data_Raw/Occurrences/calflora_data", full.names = TRUE)
calflora_raw <- map_dfr(calflora_files, ~read_csv(.x, col_types = cols(.default = "c")))

# 3) Save formatted
write.csv(calflora_raw, "Data_Raw/Formatted/calflora_formatted.csv", row.names = F)

# Offload
rm(calflora_files, calflora_raw)
gc()

##### CCH2 ##### 
### Approach
# Download occurrence from https://www.cch2.org/ using geoson files generated below, georeferenced only, exclude cultivated/captive
# There is a 1 million record download limit, so you have to do it in batches
# Code to generate the 5 GeoJSON shapefiles to download CCH2 data in batches, including 100 km buffer

# Set package/plotting options
options(tigris_use_cache = TRUE)
sf_use_s2(TRUE)  # use spherical geometry for intersections

# 1) Get California boundary 
ca_ll <- states(year = 2023, cb = TRUE, progress_bar = FALSE) |>
  st_as_sf() |>
  filter(STUSPS == "CA") |>
  st_make_valid() |>
  st_geometry()

# 2) Reproject to California Albers (meters, preserves distance) 
ca_3310 <- st_transform(ca_ll, 3310)

# 3) Create 100 km buffer 
ca_buffer <- st_buffer(ca_3310, dist = 100000) |> st_make_valid()

# 4) Simplify geometry for manageable file size
ca_simplified <- ms_simplify(ca_buffer, keep = 0.2, keep_shapes = TRUE)

# 5) Transform to WGS84 (lon/lat) for GeoJSON
ca_wgs84 <- st_transform(ca_simplified, 4326)

# 6) Split into 5 roughly equal latitude bands 
bbox <- st_bbox(ca_wgs84)
lat_breaks <- seq(bbox["ymin"], bbox["ymax"], length.out = 6)  # 5 chunks = 6 breaks

# Create horizontal slice polygons
slices <- lapply(seq_along(lat_breaks[-1]), function(i) {
  slice_poly <- st_polygon(list(rbind(
    c(bbox["xmin"], lat_breaks[i]),
    c(bbox["xmax"], lat_breaks[i]),
    c(bbox["xmax"], lat_breaks[i+1]),
    c(bbox["xmin"], lat_breaks[i+1]),
    c(bbox["xmin"], lat_breaks[i])
  ))) |> st_sfc(crs = 4326)
  
  part <- suppressWarnings(st_intersection(ca_wgs84, slice_poly))
  if (length(part) > 0 && !is.null(part)) part else NULL
})

# Clean out any empty/null slices
slices <- slices[!sapply(slices, is.null)]

# 7) Write GeoJSON files
names(slices) <- c("north", "north_central", "central", "south_central", "south")

for (name in names(slices)) {
  fname <- paste0("california_plus_100km_", name, ".geojson")
  st_write(slices[[name]], fname, driver = "GeoJSON", delete_dsn = TRUE)
  message("Saved ", fname)
}

# Optional: quick visual check 
plot(st_geometry(ca_wgs84), col = "grey90", border = "grey40")
lapply(seq_along(slices), function(i)
  plot(st_geometry(slices[[i]]), add = TRUE, border = i + 1, lwd = 2))

# Now compile the data

# Define the main folder path
cch2_folder <- "Data_Raw/Occurrences/cch2"

# Find all "occurrences.csv" files
file_paths <- list.files(
  cch2_folder,
  pattern = "occurrences\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# Read and combine (catch-all safe version)
cch2_list <- lapply(file_paths, function(f) {
  message("Reading: ", f)
  suppressWarnings(
    readr::read_csv(
      f,
      col_types = readr::cols(.default = "c"),  # force all columns as character
      progress = FALSE
    )
  )
})

# Combine safely
cch2 <- dplyr::bind_rows(cch2_list)

## I think every row of data in the CCH2 datasets are given a unique ID, so even if there were duplicates between the 3 datasets (which there are), you first need to get rid of that unique ID column to remove them. Skipping for now

# Save formatted
write.csv(cch2, "Data_Raw/Formatted/cch2_formatted.csv", row.names = F)

# Offload
rm(cch2, cch2_list)
gc()

##### CCH2 Phenology ##### 
### Approach
# Download occurrence from https://www.cch2.org/ using the geoson files generated below, georeferenced only, exclude cultivated/captive
# AND in Trait Critera, select "Open Flower" = Present

# Set package/plotting options
options(tigris_use_cache = TRUE)
sf_use_s2(TRUE)  # use spherical geometry for intersections

# Parameters 
n_slices <- 6  # <--- Change this to any number of desired latitude slices

# 1) Get California boundary 
ca_ll <- states(year = 2023, cb = TRUE, progress_bar = FALSE) |>
  st_as_sf() |>
  filter(STUSPS == "CA") |>
  st_make_valid() |>
  st_geometry()

# 2) Reproject to California Albers (meters, preserves distance)
ca_3310 <- st_transform(ca_ll, 3310)

# 3) Simplify geometry for manageable file size
ca_simplified <- ms_simplify(ca_3310, keep = 0.2, keep_shapes = TRUE)

# 4) Transform to WGS84 (lon/lat) for GeoJSON
ca_wgs84 <- st_transform(ca_simplified, 4326)

# 5) Split into n_slices equal latitude sections
bbox <- st_bbox(ca_wgs84)
lat_breaks <- seq(bbox["ymin"], bbox["ymax"], length.out = n_slices + 1)

slices <- lapply(seq_along(lat_breaks[-1]), function(i) {
  slice_poly <- st_polygon(list(rbind(
    c(bbox["xmin"], lat_breaks[i]),
    c(bbox["xmax"], lat_breaks[i]),
    c(bbox["xmax"], lat_breaks[i+1]),
    c(bbox["xmin"], lat_breaks[i+1]),
    c(bbox["xmin"], lat_breaks[i])
  ))) |> st_sfc(crs = 4326)
  
  part <- suppressWarnings(st_intersection(ca_wgs84, slice_poly))
  if (length(part) > 0 && !is.null(part)) part else NULL
})

# 6) Clean and name slices
slices <- slices[!sapply(slices, is.null)]
names(slices) <- paste0("slice_", seq_along(slices))

# 7) Write GeoJSON files 
for (name in names(slices)) {
  fname <- paste0("california_", name, ".geojson")
  st_write(slices[[name]], fname, driver = "GeoJSON", delete_dsn = TRUE)
  message("Saved ", fname)
}

# Optional: visual check 
plot(st_geometry(ca_wgs84), col = "grey90", border = "grey40")
lapply(seq_along(slices), function(i)
  plot(st_geometry(slices[[i]]), add = TRUE, border = i + 1, lwd = 2))

# Read back in
# Define the main folder path
cch2_pheno_folder <- "Data_Raw/Occurrences/cch2_pheno"

# Find all "occurrences.csv" files
file_paths <- list.files(
  cch2_pheno_folder,
  pattern = "occurrences\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# Read and combine (catch-all safe version)
cch2_pheno_list <- lapply(file_paths, function(f) {
  message("Reading: ", f)
  suppressWarnings(
    readr::read_csv(
      f,
      col_types = readr::cols(.default = "c"),  # force all columns as character
      progress = FALSE
    )
  )
})

# Combine safely
cch2_pheno <- dplyr::bind_rows(cch2_pheno_list)

#Save formatted
write.csv(cch2_pheno, "Data_Raw/Formatted/cch2_pheno_formatted.csv", row.names = F)

#Offload
rm(cch2_pheno, cch2_pheno_list)
gc()

#### Pollinators ####

##### GBIF #####

### Manual Approach
#download from https://www.gbif.org/occurrence/search 
#called the zip folder "gbif_manual.zip"
#Last downloaded Nov 20, 2024

#Have to use data.table to read it in, read.csv throws errors
gbif_raw1 = data.table::fread("Data_Raw/Occurrences/gbif/0024811-241107131044228.csv", quote="")

gbif_raw1 <- data.table::fread(
  "Data_Raw/Occurrences/gbif/0024811-241107131044228.csv",
  quote = "",  # Since you mentioned this helps
  fill = TRUE, # Handle rows with unequal length
  header = TRUE,
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)

### Automated Approach
#Supply GBIF user credentials
gbif_creds <- list(
  user = "ccosm001",
  pwd = "Wyl2smb?",
  email = "ccosm001@ucr.edu"
)

# make list of GBIF taxon keys for each desired pollinator group
name_backbone(name = "Lepidoptera") #797
name_backbone(name = "Syrphidae") #6920
name_backbone(name = "Trochilidae")#5289
# Get all bee family keys
name_backbone(name = "Apidae")
name_backbone(name = "Megachilidae")
name_backbone(name = "Halictidae")
name_backbone(name = "Andrenidae")
name_backbone(name = "Colletidae")
name_backbone(name = "Melittidae")
name_backbone(name = "Stenotritidae")
# Complete list of taxon keys for all pollinator groups
taxon_keys <- c(
  797,   # Lepidoptera (butterflies and moths)
  6920,  # Syrphidae (hover flies)
  5289,  # Trochilidae (hummingbirds)
  4334,  # Apidae (bee families here and below)
  7911,  # Megachilidae
  7908,  # Halictidae
  7901,  # Andrenidae
  7905,  # Colletidae
  4345,  # Melittidae
  7916   # Stenotritidae
)

#Generate WKT (geometry) for CA + 100km buffer. rgbif needs WKT format
# Get map sf object for CA
california_sf <- states(cb = TRUE, resolution = "500k") %>%
  filter(STUSPS == "CA")

#fix invalid geometries
california_sf <- st_make_valid(california_sf)

# Project to CA Albers for accurate 100km buffer
california_proj <- st_transform(california_sf, crs = 3310)
california_buffer <- st_buffer(california_proj, dist = 100000)

# Ensure correct winding order (counter-clockwise)
california_buffer_wgs84 <- st_reverse(california_buffer_wgs84)

#Convert to wkt
california_wkt <- st_as_text(st_geometry(california_buffer_wgs84))

# Initiate a download request in California + 100km buffer
# Last Downloaded Nov 19, 2024
key <- occ_download(
  pred_in("taxonKey",taxon_keys),  # Filter for desired groups
  pred("hasCoordinate", TRUE),        # Include only records with coordinates
  pred("occurrenceStatus", "PRESENT"), # Only presence records
  pred("geometry", california_wkt),          # Filter to California + 100km
  format = "DWCA",     # Request full darwin core
  user = gbif_creds$user,
  pwd = gbif_creds$pwd,
  email = gbif_creds$email               
)

# See status of download
occ_download_wait(key)

# Retrieve the download file and import to R (takes a few minutes; saves a zipped copy in the indicated folder path)
gbif_raw <- occ_download_get(key, path = "Data_Raw/Occurrences/gbif") %>%
  occ_download_import()

# # Or unzip the file manually and read in the occurrence.txt file
# gbif_raw = fread('Data_Raw/Occurrences/gbif/occurrence.txt', quote = "")

# Fix common issue with GBIF data
gbif_raw$scientificName <- ifelse(grepl("^BOLD", gbif_raw$scientificName), gbif_raw$species, gbif_raw$scientificName)

# Save formatted
write.csv(gbif_raw, "Data_Raw/Formatted/gbif_formatted.csv", row.names = F)

# Offload
rm(gbif_raw)
gc()

##### iDigBio #####
# Get California + 100km buffer
california_sf <- states(cb = TRUE, resolution = "500k") %>%
  filter(STUSPS == "CA") %>%
  st_make_valid()

california_proj <- st_transform(california_sf, crs = 3310)
california_buffer <- st_buffer(california_proj, dist = 100000)
california_buffer_wgs84 <- st_transform(california_buffer, crs = 4326)

# Get all counties in CA and neighboring states
all_states <- c("CA", "OR", "NV", "AZ")
all_counties <- counties(state = all_states, cb = TRUE) %>%
  st_transform(4326) %>%
  st_make_valid()

# Find counties that intersect with the 100km buffer
counties_in_buffer <- all_counties %>%
  st_filter(california_buffer_wgs84) %>%
  st_drop_geometry() %>%
  mutate(
    state_name = tolower(STATE_NAME),
    county_name = tolower(NAME)
  ) %>%
  select(state_name, county_name) %>%
  distinct() %>%
  arrange(state_name, county_name)

cat("Found", nrow(counties_in_buffer), "counties within CA + 100km buffer\n")
cat("By state:\n")
print(table(counties_in_buffer$state_name))

# Define taxa to download
taxa_queries <- list(
  Lepidoptera = list(order = "lepidoptera"),
  Syrphidae = list(family = "syrphidae"),
  Trochilidae = list(family = "trochilidae"),
  Apidae = list(family = "apidae"),
  Megachilidae = list(family = "megachilidae"),
  Halictidae = list(family = "halictidae"),
  Andrenidae = list(family = "andrenidae"),
  Colletidae = list(family = "colletidae"),
  Melittidae = list(family = "melittidae")
)

# Function to retry API calls with exponential backoff
retry_api_call <- function(func, max_attempts = 3, initial_wait = 2) {
  attempt <- 1
  while(attempt <= max_attempts) {
    result <- tryCatch({
      func()
    }, error = function(e) {
      if(grepl("500|Internal Server Error", e$message) && attempt < max_attempts) {
        wait_time <- initial_wait * (2 ^ (attempt - 1))
        cat(sprintf(" - ERROR (attempt %d/%d): %s\n", attempt, max_attempts, e$message))
        cat(sprintf("   Retrying in %d seconds...", wait_time))
        Sys.sleep(wait_time)
        return(NULL)
      } else {
        stop(e)
      }
    })
    
    if(!is.null(result)) {
      return(result)
    }
    
    attempt <- attempt + 1
  }
  
  stop("Maximum retry attempts reached")
}

# Recursive function to download records, splitting by year if needed
download_with_splitting <- function(rq, fields, date_range = NULL, indent = "    ", depth = 0, max_depth = 5) {
  
  # Add date filter if provided
  if(!is.null(date_range)) {
    date_filter <- list(type = "range")
    if(!is.null(date_range$gte)) date_filter$gte <- date_range$gte
    if(!is.null(date_range$lte)) date_filter$lte <- date_range$lte
    rq <- c(rq, list(datecollected = date_filter))
  }
  
  # Get count
  count <- retry_api_call(function() idig_count_records(rq = rq))
  
  # Print info about this chunk
  if(!is.null(date_range)) {
    label <- paste0(
      if(!is.null(date_range$gte)) paste0(">=", substr(date_range$gte, 1, 4)) else "pre",
      " to ",
      if(!is.null(date_range$lte)) paste0("<", substr(date_range$lte, 1, 4)) else "present"
    )
    cat(sprintf("%s%s: %d records", indent, label, count))
  } else {
    cat(count, "records")
  }
  
  # If count is 0, return NULL
  if(count == 0) {
    cat("\n")
    return(NULL)
  }
  
  # If count is within limit, download
  if(count <= 100000) {
    records <- retry_api_call(function() {
      idig_search_records(rq = rq, fields = fields, limit = 100000)
    })
    cat(" ✓\n")
    return(records)
  }
  
  # If we've hit max depth, warn and skip
  if(depth >= max_depth) {
    cat(" - MAX DEPTH REACHED, skipping\n")
    return(NULL)
  }
  
  # Count is > 100000, need to split
  cat(" - splitting\n")
  
  # Determine year range to split
  if(is.null(date_range) || is.null(date_range$gte)) {
    start_year <- 1800  # Default start
  } else {
    start_year <- as.integer(substr(date_range$gte, 1, 4))
  }
  
  if(is.null(date_range) || is.null(date_range$lte)) {
    end_year <- as.integer(format(Sys.Date(), "%Y"))
  } else {
    end_year <- as.integer(substr(date_range$lte, 1, 4))
  }
  
  # Calculate midpoint year
  mid_year <- floor((start_year + end_year) / 2)
  
  # Create two chunks
  chunk1_range <- list(
    gte = if(is.null(date_range) || is.null(date_range$gte)) NULL else date_range$gte,
    lte = sprintf("%d-01-01", mid_year)
  )
  
  chunk2_range <- list(
    gte = sprintf("%d-01-01", mid_year),
    lte = if(is.null(date_range) || is.null(date_range$lte)) NULL else date_range$lte
  )
  
  # Recursively download both chunks
  records_list <- list()
  
  records1 <- download_with_splitting(rq, fields, chunk1_range, 
                                      paste0(indent, "  "), depth + 1, max_depth)
  if(!is.null(records1)) {
    records_list[[length(records_list) + 1]] <- records1
  }
  
  Sys.sleep(0.5)
  
  records2 <- download_with_splitting(rq, fields, chunk2_range, 
                                      paste0(indent, "  "), depth + 1, max_depth)
  if(!is.null(records2)) {
    records_list[[length(records_list) + 1]] <- records2
  }
  
  if(length(records_list) > 0) {
    return(bind_rows(records_list))
  } else {
    return(NULL)
  }
}

# Initialize list to store results
all_records <- list()

# Download records for each taxon
for(taxon_name in names(taxa_queries)) {
  
  cat("\n=== Downloading", taxon_name, "===\n")
  
  taxon_records <- list()
  
  # Download by county to avoid API limits
  for(i in 1:nrow(counties_in_buffer)) {
    state <- counties_in_buffer$state_name[i]
    county <- counties_in_buffer$county_name[i]
    
    cat(sprintf("  %s County, %s...", tools::toTitleCase(county), toupper(state)))
    
    # Build query with county filter
    rq <- c(
      taxa_queries[[taxon_name]],
      list(
        stateprovince = state,
        county = county
      )
    )
    
    fields <- c("uuid", "scientificname", "taxonrank", 
                "order", "family", "genus", "specificepithet", "infraspecificepithet",
                "country", "stateprovince", "county", "locality",
                "geopoint", "coordinateuncertainty", 
                "datecollected", "collector", "catalognumber",
                "institutioncode", "collectioncode", "basisofrecord",
                "occurrenceid")
    
    tryCatch({
      records <- download_with_splitting(rq, fields)
      if(!is.null(records)) {
        taxon_records[[length(taxon_records) + 1]] <- records
      }
    }, error = function(e) {
      cat("  FINAL ERROR:", e$message, "\n")
    })
    
    Sys.sleep(0.5)
  }
  
  if(length(taxon_records) > 0) {
    all_records[[taxon_name]] <- bind_rows(taxon_records)
    cat("Total downloaded for", taxon_name, ":", nrow(all_records[[taxon_name]]), "records\n")
  } else {
    cat("No records downloaded for", taxon_name, "\n")
  }
  
  Sys.sleep(2)
}

# Combine all taxa into one dataframe
cat("\n=== Combining all records ===\n")
idigbio_raw <- bind_rows(all_records, .id = "query_taxon")

# Remove duplicates
idigbio_raw <- idigbio_raw %>% distinct(uuid, .keep_all = TRUE)
cat("Total unique records:", nrow(idigbio_raw), "\n")

# Save raw
write.csv(idigbio_raw, "Data_Raw/Occurrences/idigbio/idigbio_raw.csv", row.names = FALSE)

# Read in 
idigbio_raw = read.csv("Data_Raw/Occurrences/idigbio/idigbio_raw.csv")

# Save formatted
write.csv(idigbio_raw, "Data_Raw/Formatted/idigbio_formatted.csv", row.names = F)

# Offload
rm(idigbio_raw)
gc()

##### Ecdysis (Formerly SCAN) #####
### Approach
# Scan is now on ecdysis: https://ecdysis.org/collections/search/index.php
# Same download structure as CCH2
# Create a geoson shape of ca + buffer for ecdysis, copy and paste in the polygon search field

options(tigris_use_cache = TRUE)
sf_use_s2(TRUE)  # use spherical geometry for intersections

# 1) Get California boundary 
ca_ll <- states(year = 2023, cb = TRUE, progress_bar = FALSE) |>
  st_as_sf() |>
  filter(STUSPS == "CA") |>
  st_make_valid() |>
  st_geometry()

# 2) Reproject to California Albers (meters, preserves distance)
ca_3310 <- st_transform(ca_ll, 3310)

# 3) Create 100 km buffer 
ca_buffer <- st_buffer(ca_3310, dist = 100000) |> 
  st_make_valid()

# 4) Simplify geometry for manageable file size 
ca_simplified <- ms_simplify(ca_buffer, keep = 0.2, keep_shapes = TRUE)

# 5) Transform to WGS84 (lon/lat) for GeoJSON
ca_wgs84 <- st_transform(ca_simplified, 4326)

# 6) Write single GeoJSON file
fname <- "california_plus_100km.geojson"
st_write(ca_wgs84, fname, driver = "GeoJSON", delete_dsn = TRUE)
message("Saved ", fname)

# Optional: quick visual check 
plot(st_geometry(ca_wgs84), col = "grey90", border = "grey40", main = "California + 100km Buffer")

# Read in 
ecdysis_raw = fread('/Users/chriscosma/Desktop/Morpho Code Refinements 10:2025/input/Data_Raw/Occurrences/ecdysis/SymbOutput_2025-10-17_175217_DwC-A/occurrences.csv')

# Save formatted
write.csv(ecdysis_raw, "Data_Raw/Formatted/ecdysis_formatted.csv", row.names = F)

# Offload
rm(ecdysis_raw)
gc()

##### MPG #####

#Records supplied by Steve Nanz of Moth Photographers' Group on Nov 13, 2024
mpg_raw = read.csv("Data_Raw/Occurrences/mpg/mpg_raw.csv")

# Save formatted
write.csv(mpg_raw, "Data_Raw/Formatted/mpg_formatted.csv", row.names = F)

# Offload
rm(mpg_raw)
gc()

##### Chesshire #####

# Download from https://figshare.com/projects/Completeness_analyses_for_over_3000_United_States_bee_species_identifies_persistent_data_gaps/138673 

chesshire_raw <- read.csv("Data_Raw/Occurrences/chesshire/contiguousRecords_high_Only.csv") 

# Save formatted
write.csv(chesshire_raw, "Data_Raw/Formatted/chesshire_formatted.csv", row.names = F)

# Offload
rm(chesshire_raw)
gc()
