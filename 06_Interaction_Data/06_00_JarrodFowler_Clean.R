#### Load Packages ####
library(tidyverse)
library(readxl)

#### Load Data ####
jarrodFowler_specialistBees <- read_excel("Data_Raw/Interactions/JarrodFowler/Fowler_SpecialistBees.xlsx")

#### Clean Table ####
pollinator_col <- 2 
host_plant_col <- 18

# Extract relevant columns
pollinator_plant_data <- jarrodFowler_specialistBees %>%
  # Select only pollinator species and host plant columns
  dplyr::select(pollinator = all_of(pollinator_col), 
         host_plants = all_of(host_plant_col)) %>%
  # Remove rows with missing data
  filter(!is.na(pollinator) & !is.na(host_plants)) %>%
  # Remove empty strings
  filter(pollinator != "" & host_plants != "")

# Function to split multiple plant species and create individual rows
# handles various separators that are  used (comma, semicolon, "and", etc.)
expand_plant_interactions <- function(df) {
  df %>%
    # Separate host_plants by common delimiters
    separate_rows(host_plants, sep = ",|;|:|\\band\\b") %>%
    # Clean up whitespace
    mutate(host_plants = str_trim(host_plants)) %>%
    # Remove any remaining empty values
    filter(host_plants != "" & !is.na(host_plants))
}

# Create the final dataframe
interaction_df <- pollinator_plant_data %>%
  expand_plant_interactions() %>%
  # Rename columns to match desired format
  rename(higher = pollinator,
         lower = host_plants) %>%
  # Remove subgenus names (text within parentheses) from pollinator names
  mutate(higher = str_remove(higher, "\\s*\\([^)]+\\)")) %>%
  # Add interaction type
  mutate(interactionTypeName = "collectsPollenOf") %>%
  # Select columns in desired order
  dplyr::select(higher, lower, interactionTypeName) %>%
  # Remove duplicate rows (same pollinator-plant combination)
  distinct()

# Save the result to CSV
write.csv(interaction_df, "Data_Raw/Interactions/JarrodFowler/jarrodFowler_scraped_interactions.csv", row.names = FALSE)
