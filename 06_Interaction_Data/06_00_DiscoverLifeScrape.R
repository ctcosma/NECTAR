
#### Load Packages ####
library(rvest)
library(stringr)
library(dplyr)
library(purrr)

#### Functions ####
# Function for scraping interaction record IDs from Discover Life
get_record_ids <- function(species_name) {
  # Format URL
  url <- paste0("https://www.discoverlife.org/mp/20q?search=", species_name)
  
  # add delay to prevent server overload
  Sys.sleep(0.2)
  
  tryCatch({
    # Read and parse the HTML
    page <- read_html(url)
    
    # Get all rows from the Hosts table that contain record links
    host_rows <- page %>%
      html_nodes("tr") %>%
      # Keep only rows that don't contain "Withheld" as the scientific name
      .[!str_detect(html_text(.), "Withheld")]
    
    # Extract links from these rows
    links <- host_rows %>%
      html_nodes("a") %>%
      html_attr("href") %>%
      # Keep only links that match the pattern for record IDs
      .[str_detect(., "/mp/20l\\?id=")]
    
    # Extract IDs from the links
    record_ids <- links %>%
      # Split on 'id=' and take the second part
      str_extract("(?<=id=).*") %>%
      # Split ids that are separated by semicolons
      str_split(";") %>%
      unlist() %>%
      unique()
    
    return(record_ids)
    
  }, error = function(e) {
    warning(paste("Error fetching", url, ":", e$message))
    return(character(0))
  })
}

# Function to scrape a single record
get_record_data <- function(record_id) {
  # what url are we querying
  url <- paste0("https://www.discoverlife.org/mp/20p?see=", record_id)
  
  # Add delay
  Sys.sleep(0.2)
  
  tryCatch({
    # Read the page
    page <- read_html(url)
    
    tables <- html_nodes(page, "table")
    
    # Initialize record data
    record_data <- list()
    
    # Process all tables
    for(table in tables) {
      data_rows <- html_nodes(table, "tr")
      
      for(row in data_rows) {
        cells <- html_nodes(row, "td")
        if(length(cells) >= 2) {
          field <- cells[1] %>% html_text() %>% trimws()
          value <- cells[2] %>% html_text() %>% trimws()
          
          # Only add if we have a valid field name and it's not already in record_data
          if(nchar(field) > 0 && !field %in% c("", "MAP") && 
             !(field %in% names(record_data))) {
            record_data[[field]] <- value
          }
        }
      }
    }
    
    # Add record ID
    record_data$record_id <- id
    records_data[[id]] <- record_data
    
  }, error = function(e) {
    warning(paste("Error processing record", id, ":", e$message))
    records_data[[id]] <- list(record_id = id, error = e$message)
  })
}

#### Load Data ####
## read in checklists
# Plant checklist
plant_checklist <- read.csv("Data_Clean/Species_Checklists/plant_checklist.csv")  %>%
  select(scientificName) %>%
  mutate(genera = str_extract(scientificName, "[A-Z][a-z]*")) %>%
  mutate(scientificName = str_extract(scientificName, "^[A-Z][a-z]*(?:\\s[a-z]+)?"))

# Pollinator checklists
bee_checklist <- read.csv("Data_Clean/Species_Checklists/bee_checklist.csv")
hoverfly_checklist <- read.csv("Data_Clean/Species_Checklists/hoverfly_checklist.csv")
butterfly_checklist <- read.csv("Data_Clean/Species_Checklists/butterfly_checklist.csv")
moth_checklist <- read.csv("Data_Clean/Species_Checklists/moth_checklist.csv")

higher_checklists <- bind_rows(butterfly_checklist, moth_checklist, bee_checklist, hoverfly_checklist) %>%
  select(scientificName) %>%
  mutate(genera = str_extract(scientificName, "[A-Z][a-z]*")) %>%
  mutate(scientificName = str_extract(scientificName, "^[A-Z][a-z]*(?:\\s[a-z]+)?"))

#### Discover Life Parse/Scrape ####
### Get record IDs for interaction records
# collapse species names, prepping for web scraping
parsing_names <- gsub(" ", "+", higher_checklists$scientificName)

# Initialize empty list for storing record IDs
record_IDs <- list()

# Process each species and store in the list
for(sp in parsing_names) {
  cat("Processing", sp, "\n")  # Progress message
  record_IDs[[sp]] <- get_record_ids(sp)
}

# remove NA from all items of list
record_IDs_clean <- map(record_IDs, ~.x[!is.na(.x)])

saveRDS(record_IDs_clean, "Data_Raw/Interactions/DiscoverLife/recordIDs.rds")

### Now extract data from table of interaction records
for (sp in names(record_IDs_clean)){
  # progress indicator
  cat( "\n", "Extracting", length(record_IDs_clean[[sp]]), "records for", sp, "\n")
  
  # Process each record, then save species 
  records_data <- list()
  
  for(id in record_IDs_clean[[sp]]) {
    cat(".", sep = "")  # Progress indicator
    record_data <- get_record_data(id)
    records_data[[id]] <- record_data
  }
  
  # Save species records
  saveRDS(records_data, 
          file = file.path("Data_Raw/Interactions/DiscoverLife/species_records_raw", paste0(sp, ".rds")))
}


