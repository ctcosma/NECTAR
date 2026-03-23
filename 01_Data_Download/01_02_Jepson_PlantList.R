#### Load packages ####
library(rvest)
library(dplyr)
library(stringr)

#### Functions ####
# Scrape data from a single index page on Jepson eFlora
scrape_jepson_index <- function(letter) {
  url <- paste0("https://ucjeps.berkeley.edu/eflora/eflora_index.php?index=", letter)
  
  cat("Scraping index:", letter, "\n")
  
  # Read the HTML
  page <- read_html(url)
  
  # Find the table within the eFloraTable div
  table <- page %>%
    html_element("div.eFloraTable table")
  
  # Extract all rows (skip header row)
  rows <- table %>%
    html_elements("tr") %>%
    tail(-1)  # Remove header row
  
  # Extract data from each row
  data_list <- lapply(rows, function(row) {
    cells <- row %>% html_elements("td")
    
    if (length(cells) == 3) {
      # Extract scientific name (clean up HTML tags and special characters)
      sci_name <- cells[1] %>%
        html_text2() %>%
        str_replace_all("\\s+", " ") %>%
        str_trim()
      
      # Extract common name
      common_name <- cells[2] %>%
        html_text2() %>%
        str_trim()
      
      # Extract status
      status <- cells[3] %>%
        html_text2() %>%
        str_trim()
      
      return(data.frame(
        scientific_name = sci_name,
        common_name = common_name,
        status = status,
        index_letter = letter,
        stringsAsFactors = FALSE
      ))
    } else {
      return(NULL)
    }
  })
  
  # Combine into dataframe
  do.call(rbind, data_list)
}

#### Scrape Jepson web data ####
letters_to_scrape <- LETTERS

# Scrape all pages with error handling
all_data <- lapply(letters_to_scrape, function(letter) {
  tryCatch({
    Sys.sleep(1)  # Be polite - wait 1 second between requests
    scrape_jepson_index(letter)
  }, error = function(e) {
    cat("Error with letter", letter, ":", e$message, "\n")
    return(NULL)
  })
})

# Combine all dfs
jepson_plants <- do.call(rbind, all_data)

# Clean up the data
jepson_plants <- jepson_plants %>%
  mutate(
    # Extract accepted name from synonyms (text within "Under ...")
    accepted_name = str_extract(scientific_name, "(?<=\\(Under )([^)]+)(?=\\))"),
    
    # Check if entry has a key
    has_key = str_detect(scientific_name, "Key to"),
    
    # Clean scientific name by removing "↳ Key to ..." and "(Under ...)" text
    scientific_name_clean = scientific_name %>%
      str_remove("↳\\s*Key to.*$") %>%  # Remove key notation
      str_remove("\\(Under [^)]+\\)") %>%  # Remove synonym reference
      str_trim(),
    
    # Replace empty strings with NA
    common_name = ifelse(common_name == "", NA, common_name),
    status = ifelse(status == "", NA, status),
    accepted_name = ifelse(is.na(accepted_name) | accepted_name == "", NA, accepted_name)
  ) %>%
  # Reorder columns for clarity
  select(scientific_name_clean, scientific_name, common_name, status, 
         accepted_name, has_key, index_letter)

jepson_plants_cleaned <- jepson_plants %>%
  select(scientific_name_clean, common_name, status, accepted_name) %>%
  rename(scientificName = scientific_name_clean,
         acceptedName = accepted_name,
         commonName = common_name)

# Clean up names for matching with parsed names
jepson_plants_cleaned <- jepson_plants_cleaned %>%
  mutate(
    scientificNameClean = stringr::str_remove_all(scientificName, paste(c("subsp. ", "var. ", "f. "), collapse = "|")),
    acceptedNameClean = stringr::str_remove_all(acceptedName, paste(c("subsp. ", "var. ", "f. "), collapse = "|"))
  )

#### Output ####
# Save to CSV
write.csv(jepson_plants_cleaned, "Data_Raw/Species_Traits_Attributes/jepson_eflora_plants.csv", row.names = FALSE)

