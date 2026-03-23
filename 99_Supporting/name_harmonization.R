#### Summary ####
## Core harmonization pipeline
# names_harmonize:          Send names to gnverifier
# names_choose:             Select best match from gnverifier output
# clean_harmonized_names:   Parse and clean gnverifier output names
# correct_names:            Apply manual name corrections
# harmonize_names:          Wrapper for the full harmonization pipeline
#
## Plant-specific
# apply_jepson_taxonomy:    Update plant names using Jepson eFlora
#
## Taxonomy
# update_taxonomy:          Extract higher taxonomy from ClassificationPath

#### Required packages ####
library(dplyr)
library(data.table)
library(stringr)
library(bdc)

#### Core Harmonization Functions ####
# Send a list of names to gnverifier and return full match results
names_harmonize <- function(names, 
                            names_file = "Temp/names.tsv", 
                            names_harm_file = "Temp/names_harmonized.tsv",
                            gnv_path = dirname(Sys.which("gnverifier"))) {
  
  # Save names as tsv
  write.table(names, names_file, sep = "\t", row.names = FALSE, col.names = FALSE, 
              quote = FALSE, fileEncoding = "UTF-8")
  
  # Run gnverifier (https://github.com/gnames/gnverifier)
  system(paste(file.path(gnv_path, "gnverifier"), "-M", shQuote(names_file), ">", shQuote(names_harm_file)))
  
  # Read in results
  names_harm <- fread(names_harm_file, sep = ",", header = TRUE, encoding = "UTF-8")
  
  return(names_harm)
}

# Select best match from gnverifier output based on higher taxonomy and sort score
names_choose <- function(names_harm, higher_tax = "Plantae|Lepidoptera|Hymenoptera|Diptera|Trochilidae") {
  
  names_harm_withHigher <- names_harm %>%
    filter(grepl(higher_tax, ClassificationPath)) 
  
  names_harm_select <- names_harm_withHigher %>% 
    group_by(ScientificName) %>% 
    arrange(SortScore) %>%
    {
      bind_rows(
        filter(., Kind == "BestMatch") %>% slice_head(n = 1),
        anti_join(., filter(., Kind == "BestMatch"), by = "ScientificName") %>% slice_head(n = 1)
      )
    } %>%
    ungroup()
  
  return(names_harm_select)
}

# Parse and clean CurrentName from gnverifier output using bdc_clean_names
clean_harmonized_names <- function(names_harm_select) {
  
  # Re-parse CurrentName (has authorship attached) to get clean names
  parse_names <-
    bdc_clean_names(sci_names = names_harm_select$CurrentName, save_outputs = FALSE)
  
  names_harm_clean <- names_harm_select %>%
    mutate(names_harm = parse_names$names_clean) %>%
    mutate(names_harm = ifelse(is.na(names_harm), MatchedCanonical, names_harm)) %>%
    mutate(names_harm = ifelse(MatchType == "NoMatch", NA, names_harm)) %>%
    mutate(names_harm = stringr::str_remove_all(names_harm, paste(c("subsp. ", "var. ", "f. "), collapse = "|"))) 
  
  return(names_harm_clean)
}

# Apply manual name corrections from a CSV file (actions: "change" or "remove" as listed in file)
correct_names <- function(dataset, file = "Data_Raw/taxonomy_corrections.csv", column) {
  
  corrections <- read.csv(file, stringsAsFactors = FALSE)
  col_name <- deparse(substitute(column))
  
  # Apply "change" corrections
  changes <- corrections[corrections$action == "change", ]
  if (nrow(changes) > 0) {
    for (i in 1:nrow(changes)) {
      dataset[[col_name]][dataset[[col_name]] == changes$old_name[i]] <- changes$new_name[i]
    }
  }
  
  # Apply "remove" corrections
  removals <- corrections[corrections$action == "remove", ]
  if (nrow(removals) > 0) {
    dataset <- dataset[!dataset[[col_name]] %in% removals$old_name, ]
  }
  
  return(dataset)
}

# Wrapper: run full harmonization pipeline for one taxonomic group
harmonize_names <- function(names,
                            higher_tax,
                            names_file,
                            names_harm_file,
                            gnv_path = "/home/tb625/bin/",
                            taxonomy_corrections_file = "Data_Raw/taxonomy_corrections.csv") {
  
  names_harm <- names_harmonize(names,
                                names_file = names_file,
                                names_harm_file = names_harm_file,
                                gnv_path = gnv_path)
  
  names_harm_select <- names_choose(names_harm, higher_tax = higher_tax)
  
  names_harm <- clean_harmonized_names(names_harm_select)
  
  names_harm <- correct_names(names_harm, file = taxonomy_corrections_file, column = names_harm)
  
  return(names_harm)
}

#### Plant specific helper function ####
# Update harmonized plant names using Jepson eFlora as taxonomic authority
apply_jepson_taxonomy <- function(names_harm_plants,
                                  jepson_file = "Data_Raw/Species_Traits_Attributes/jepson_eflora_plants.csv") {
  
  # Load and prepare Jepson names
  jepson_names <- read.csv(jepson_file) %>%
    dplyr::select(scientificNameClean, acceptedNameClean, status) %>%
    rename(names_harm_jepson = acceptedNameClean) %>%
    filter(!grepl(", in part", status),
           !status %in% c("Illegitimate name",
                          "Invalid name",
                          "Misapplied name",
                          "Noted name",
                          "Unabridged misapplied name")) %>%
    mutate(names_harm_jepson = ifelse(!is.na(status) & is.na(names_harm_jepson),
                                      scientificNameClean,
                                      names_harm_jepson))
  
  # Remove synonyms that map to multiple accepted names (i.e. species that have been split)
  duplicate_synonyms <- jepson_names %>%
    filter(status == "Synonym") %>%
    group_by(scientificNameClean) %>%
    filter(n() > 1 & n_distinct(names_harm_jepson) > 1) %>%
    ungroup()
  
  jepson_names <- jepson_names %>%
    anti_join(duplicate_synonyms, by = c("scientificNameClean", "names_harm_jepson", "status"))
  
  # Four-pass matching against Jepson
  names_harm_plants <- names_harm_plants %>%
    # Pass 1: original name
    left_join(jepson_names, by = join_by("ScientificName" == "scientificNameClean")) %>%
    rename(names_harm_jepson_direct = names_harm_jepson,
           status_direct = status) %>%
    # Pass 2: original name without infraspecific epithet
    mutate(ScientificName_GSP = sapply(strsplit(ScientificName, "\\s+"), function(x) paste(x[1:2], collapse = " "))) %>%
    left_join(jepson_names, by = join_by("ScientificName_GSP" == "scientificNameClean")) %>%
    rename(names_harm_jepson_direct_GSP = names_harm_jepson,
           status_direct_GSP = status) %>%
    # Pass 3: gnverifier harmonized name
    left_join(jepson_names, by = join_by("names_harm" == "scientificNameClean")) %>%
    rename(names_harm_jepson_indirect = names_harm_jepson,
           status_indirect = status) %>%
    # Pass 4: harmonized name without infraspecific epithet
    mutate(names_harm_GSP = sapply(strsplit(names_harm, "\\s+"), function(x) paste(x[1:2], collapse = " "))) %>%
    left_join(jepson_names, by = join_by("names_harm_GSP" == "scientificNameClean")) %>%
    rename(names_harm_jepson_indirect_GSP = names_harm_jepson,
           status_indirect_GSP = status) %>%
    # Select best available match
    mutate(names_harm = coalesce(names_harm_jepson_direct, names_harm_jepson_direct_GSP,
                                 names_harm_jepson_indirect, names_harm_jepson_indirect_GSP)) %>%
    dplyr::select(-c(names_harm_jepson_direct, names_harm_jepson_direct_GSP,
                     names_harm_jepson_indirect, names_harm_jepson_indirect_GSP,
                     status_direct, status_direct_GSP, status_indirect, status_indirect_GSP,
                     ScientificName_GSP, names_harm_GSP)) %>%
    unique()
  
  return(names_harm_plants)
}

#### Higher Taxonomy ####
# Extract higher taxonomy columns from ClassificationPath column
update_taxonomy <- function(df, tax) {
  dt <- as.data.table(df)
  
  tax_lists <- strsplit(dt$ClassificationPath, "|", fixed = TRUE)
  
  if (tax == "plant") {
    dt[, `:=`(
      kingdom = sapply(tax_lists, function(x) if("Plantae" %in% x) "Plantae" else NA_character_),
      phylum = sapply(tax_lists, function(x) grep("phyta$", x, value = TRUE)[1]),
      class = sapply(tax_lists, function(x) grep("opsida$", x, value = TRUE)[1]),
      order = sapply(tax_lists, function(x) grep("ales$", x, value = TRUE)[1]),
      superfamily = sapply(tax_lists, function(x) grep("acea$", x, value = TRUE)[1]),
      family = sapply(tax_lists, function(x) grep("aceae$", x, value = TRUE)[1]),
      tribe = sapply(tax_lists, function(x) {
        matches <- grep("eae$", x, value = TRUE)
        matches <- matches[!grepl("aceae$|oideae$|acea$", matches)]
        if (length(matches) > 0) matches[1] else NA_character_
      }),
      subtribe = sapply(tax_lists, function(x) {
        matches <- grep("inae$", x, value = TRUE)
        matches <- matches[!grepl("phytina$", matches)]
        if (length(matches) > 0) matches[1] else NA_character_
      })
    )]
  } else if (tax == "animal") {
    dt[, `:=`(
      kingdom = sapply(tax_lists, function(x) if("Animalia" %in% x) "Animalia" else NA_character_),
      phylum = sapply(tax_lists, function(x) grep("arthropoda$|chordata$|mollusca$", x, value = TRUE, ignore.case = TRUE)[1]),
      class = sapply(tax_lists, function(x) grep("(insecta|arachnida|malacostraca|aves|mammalia|reptilia|amphibia)$", x, value = TRUE, ignore.case = TRUE)[1]),
      order = sapply(tax_lists, function(x) {
        matches <- grep("(ida|iformes|ptera|optera)$", x, value = TRUE, perl = TRUE)
        if (length(matches) > 0) matches[1] else NA_character_
      }),
      superfamily = sapply(tax_lists, function(x) grep("oidea$", x, value = TRUE)[1]),
      family = sapply(tax_lists, function(x) grep("idae$", x, value = TRUE)[1]),
      tribe = sapply(tax_lists, function(x) grep("ini$", x, value = TRUE)[1]),
      subtribe = sapply(tax_lists, function(x) {
        matches <- grep("ina$", x, value = TRUE)
        matches <- matches[!grepl("(inae|ini|oidea)$", matches)]
        if (length(matches) > 0) matches[1] else NA_character_
      })
    )]
  }
  
  name_parts <- strsplit(dt$names_harm, " ", fixed = TRUE)
  
  dt[, `:=`(
    genus = sapply(name_parts, function(x) if(length(x) >= 1) x[1] else NA_character_),
    specificEpithet = sapply(name_parts, function(x) if(length(x) >= 2) x[2] else NA_character_),
    infraspecificEpithet = sapply(name_parts, function(x) if(length(x) >= 3) x[3] else NA_character_)
  )]
  
  return(as.data.frame(dt))
}