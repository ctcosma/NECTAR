# Script to parse the Big Book of Hymenoptera

#### Load Packages ####
library(stringr)
library(dplyr)
library(tidyr)

#### Functions ####
parse_hymenoptera_book <- function(file_path) {
  
  # Read the file
  text <- readLines(file_path, warn = FALSE)
  
  # Initialize tracking variables
  current_genus <- NA
  current_species_group <- NA
  
  # Initialize result lists
  interactions_list <- list()
  species_list <- list()
  parasites_list <- list()
  
  # Loop through lines
  i <- 1
  while (i <= length(text)) {
    line <- text[i]
    
    # Skip empty lines
    if (nchar(trimws(line)) == 0) {
      i <- i + 1
      next
    }
    
    # Check for Genus header
    if (grepl("^Genus [A-Z]", line)) {
      current_genus <- str_extract(line, "(?<=Genus )\\w+")
      i <- i + 1
      next
    }
    
    # Check for SPECIES GROUP header
    if (grepl("^SPECIES GROUP", line)) {
      # Extract the group name (everything after "SPECIES GROUP ")
      group_match <- str_extract(line, "(?<=SPECIES GROUP )\\w+")
      if (!is.na(group_match)) {
        current_species_group <- group_match
      }
      i <- i + 1
      next
    }
    
    # Check for species entry
    # Species lines: start with lowercase, have Author name, then period
    # Example: "punctatissima Michener. Southern Ariz..."
    # Example: "fulgidus fulgidus Swenk. B. C., Mont..."
    if (grepl("^[a-z]+\\s+", line) && 
        grepl("^[a-z]+(?:\\s+[a-z]+)?\\s+[A-Z][a-z]+\\.", line) &&
        !grepl("^(Type-species:|Biology:|Taxonomy:|Revision:|Ecology:|Predator:|Parasite:)", line, ignore.case = TRUE)) {
      
      # BEFORE processing this potential species, check previous-line context:
      # If the previous line ends with a citation marker (male/female symbol codes like "2.", "'a2.", etc.)
      # or if the previous line is a Predator/Parasite line, then this lowercase+Author line is likely NOT a new Colletes species.
      prev_line_outer <- if (i > 1) trimws(text[i - 1]) else ""
      ends_with_citation_marker_outer <- grepl("(2,\\s*6\\.|2,\\s*9\\.|6\\.|2\\.|'a2\\.|'b0\\.|9\\.|d\\.|3\\.|@\\.)\\s*$", prev_line_outer)
      prev_is_predator_or_parasite <- grepl("^(Predator:|Parasite:)", prev_line_outer, ignore.case = TRUE)
      if (ends_with_citation_marker_outer || prev_is_predator_or_parasite) {
        cat("  SKIP: Previous line context indicates this is not a new Colletes species (prev_line):", prev_line_outer, "\n")
        i <- i + 1
        next
      }
      
      cat("\n=== POTENTIAL SPECIES at line", i, ":", substr(line, 1, 60), "...\n")
      
      if (is.na(current_genus)) {
        cat("  SKIP: No current genus\n")
        i <- i + 1
        next
      }
      
      # Check if this might be a continuation from previous line
      # If previous line exists and ends with a capitalized word (genus name),
      # this is likely a continuation (e.g., "Philanthus" on prev line, "albopilosus" on this line)
      if (i > 1) {
        prev_line <- text[i - 1]
        if (nchar(trimws(prev_line)) > 0 && grepl("[A-Z][a-z]+\\s*$", prev_line)) {
          cat("  SKIP: Previous line ends with capitalized word:", substr(prev_line, max(1, nchar(prev_line)-30), nchar(prev_line)), "\n")
          # Previous line ends with a capitalized word - likely a genus name
          # This line is probably a species epithet continuation, not a new entry
          i <- i + 1
          next
        }
      }
      
      # Extract species name components
      species_match <- str_match(line, "^([a-z]+)(?:\\s+([a-z]+))?\\s+([A-Z][a-z]+)\\.")
      
      if (is.na(species_match[1])) {
        cat("  SKIP: Regex match failed\n")
        i <- i + 1
        next
      }
      
      species_epithet <- species_match[2]
      subspecies_epithet <- species_match[3]  # May be NA
      initial_author <- species_match[4]
      
      # Determine full species name
      if (!is.na(subspecies_epithet)) {
        full_species <- paste(current_genus, species_epithet, subspecies_epithet)
      } else {
        full_species <- paste(current_genus, species_epithet)
      }
      
      cat("  PROCESSING:", full_species, "\n")
      
      # Accumulate all text for this species until citation line
      species_text <- line
      j <- i + 1
      citation_line <- NA
      
      # --- Inner loop: accumulate text for this species until we hit a citation block OR next species header ---
      while (j <= length(text)) {
        current_line <- text[j]
        
        # Normalize a trimmed version for pattern checks
        tline <- trimws(current_line)
        
        # --- Patterns ---
        # species header pattern used in the outer loop to detect species entries:
        # e.g., "punctatissima Michener."  OR "fulgidus fulgidus Swenk."
        species_header_pattern <- "^[a-z]+(?:\\s+[a-z]+)?\\s+[A-Z][a-z]+\\."
        
        # robust citation pattern: Genus species (or subspecies) ... , YEAR.
        citation_pattern <- paste0(
          "^[\\\\\\s]*",                    # allow leading backslash or spaces
          current_genus, "\\s+",            # genus (exact)
          "(", species_epithet,             # species epithet
          if (!is.na(subspecies_epithet)) paste0("|", subspecies_epithet) else "", ")",
          "\\b.*?(,\\s*\\d{4}\\.)"          # require comma-year-period somewhere after name
        )
        
        # fallback generic citation pattern (handles odd formatting)
        generic_citation_pattern <- paste0(
          "^[\\\\\\s]*", current_genus, "\\s+[a-z]+\\s+[A-Z][a-z]+,\\s*\\d{4}\\."
        )
        
        # line that begins with genus (uppercase) — usually a citation or synonym line
        starts_with_genus <- grepl(paste0("^[\\\\\\s]*", current_genus, "\\b"), current_line)
        
        # does this line look like a species header (implied-genus, lowercase epithet + Author + period)?
        is_species_header <- grepl(species_header_pattern, tline) &&
          !grepl("^(Type-species:|Biology:|Taxonomy:|Revision:|Ecology:|Predator:|Parasite:)", tline, ignore.case = TRUE)
        
        # is it a citation line?
        is_citation <- grepl(citation_pattern, current_line)
        if (!is_citation) {
          is_citation <- grepl(generic_citation_pattern, current_line)
        }
        
        # --- Classification decisions ---
        # 1) If this line is a citation (Colletes X Author, YEAR.) -> accumulate as citation/synonym
        if (is_citation || (starts_with_genus && grepl(",", current_line))) {
          # treat as citation/synonym line
          if (is.na(citation_line)) {
            citation_line <- current_line
            cat("  Found citation:", substr(current_line, 1, 80), "...\n")
          } else {
            cat("  Found additional citation/synonym:", substr(current_line, 1, 80), "...\n")
          }
          species_text <- paste(species_text, current_line)
          j <- j + 1
          next
        }
        
        # 2) If it starts with the genus (but didn't match the strict citation pattern),
        #    it's still most likely a synonym/citation (e.g., missing year) -> accumulate
        if (starts_with_genus) {
          cat("  Detected genus-leading line (treating as synonym/citation):", substr(current_line, 1, 80), "...\n")
          species_text <- paste(species_text, current_line)
          j <- j + 1
          next
        }
        
        # 3) If it matches the species-header pattern (lowercase epithet + Author.), treat it as the NEXT species header
        #    i.e., stop accumulating for the current species so outer loop can handle the new species.
        if (is_species_header) {
          cat("  Found next species header (stop here):", substr(current_line, 1, 80), "...\n")
          break
        }
        
        # 4) Section headers or biology/predator keywords also end the species block
        if (grepl("^(Genus|TRIBE|SUBFAMILY|Family|SPECIES GROUP|Biology:|Predator:|Parasite:)", current_line, ignore.case = TRUE)) {
          cat("  Found section boundary / keyword (stop):", substr(current_line, 1, 80), "...\n")
          break
        }
        
        # 5) Otherwise, it's continuation/descriptive text — accumulate
        if (nchar(tline) > 0) {
          species_text <- paste(species_text, current_line)
        }
        
        j <- j + 1
      } # end inner while
      
      
      # Update loop counter
      # If we found a citation, we've already moved j past it
      # If we didn't find a citation, j points to the next line to process
      i <- j
      
      # Extract author with year from citation line
      author_with_year <- NA
      if (!is.na(citation_line)) {
        # Try species epithet first
        author_pattern <- paste0(species_epithet, "\\s+([A-Z][^,]+,\\s*\\d{4})")
        author_match <- str_match(citation_line, author_pattern)
        
        # If that didn't work and we have subspecies, try subspecies
        if (is.na(author_match[1]) && !is.na(subspecies_epithet)) {
          author_pattern <- paste0(subspecies_epithet, "\\s+([A-Z][^,]+,\\s*\\d{4})")
          author_match <- str_match(citation_line, author_pattern)
        }
        
        if (!is.na(author_match[1])) {
          author_with_year <- trimws(author_match[2])
        }
      }
      
      # Extract distribution
      # Distribution is between first period and keywords (Pollen/Parasite/Ecology)
      # If no keywords, distribution goes until we see something that looks like a citation or Biology/Taxonomy line
      distribution <- NA
      first_period <- regexpr("\\.", species_text)[1]
      pollen_pos <- regexpr("Pollen:", species_text, ignore.case = TRUE)[1]
      parasite_pos <- regexpr("Parasite:", species_text, ignore.case = TRUE)[1]
      ecology_pos <- regexpr("Ecology:", species_text, ignore.case = TRUE)[1]
      biology_pos <- regexpr("Biology:", species_text, ignore.case = TRUE)[1]
      
      keyword_positions <- c(pollen_pos, parasite_pos, ecology_pos, biology_pos)
      keyword_positions <- keyword_positions[keyword_positions > 0]
      
      if (first_period > 0) {
        if (length(keyword_positions) > 0) {
          # End at first keyword
          earliest_keyword <- min(keyword_positions)
          if (earliest_keyword > first_period) {
            distribution <- substr(species_text, first_period + 1, earliest_keyword - 1)
            distribution <- trimws(distribution)
          }
        } else {
          # No keywords - distribution goes until citation line
          if (!is.na(citation_line)) {
            citation_in_text <- regexpr(paste0(current_genus, "\\s+", species_epithet, "\\s+[A-Z]"), species_text)
            if (citation_in_text < 0 && !is.na(subspecies_epithet)) {
              citation_in_text <- regexpr(paste0(current_genus, "\\s+", subspecies_epithet, "\\s+[A-Z]"), species_text)
            }
            if (citation_in_text > first_period) {
              distribution <- substr(species_text, first_period + 1, citation_in_text - 1)
              distribution <- trimws(distribution)
            }
          } else {
            # No citation found either - take everything after first period
            distribution <- substr(species_text, first_period + 1, nchar(species_text))
            distribution <- trimws(distribution)
          }
        }
      }
      
      # Extract parasite information
      parasite_info <- ""
      if (parasite_pos > 0) {
        parasite_start <- parasite_pos + 9  # Skip "Parasite:"
        
        # Find end of parasite section
        parasite_end <- nchar(species_text)
        next_keyword_positions <- c(pollen_pos, ecology_pos)
        next_keyword_positions <- next_keyword_positions[next_keyword_positions > parasite_pos]
        
        if (length(next_keyword_positions) > 0) {
          parasite_end <- min(next_keyword_positions) - 1
        } else if (!is.na(citation_line)) {
          # Try multiple patterns to find the citation
          citation_in_text <- regexpr(paste0(current_genus, "\\s+", species_epithet, "\\s+[A-Z]"), species_text)
          
          if (citation_in_text < 0 && !is.na(subspecies_epithet)) {
            citation_in_text <- regexpr(paste0(current_genus, "\\s+", subspecies_epithet, "\\s+[A-Z]"), species_text)
          }
          
          if (citation_in_text < 0) {
            citation_search <- substr(citation_line, 1, min(30, nchar(citation_line)))
            citation_in_text <- regexpr(citation_search, species_text, fixed = TRUE)
          }
          
          if (citation_in_text > 0 && citation_in_text > parasite_start) {
            parasite_end <- citation_in_text - 1
          }
        }
        
        parasite_info <- substr(species_text, parasite_start, parasite_end)
        parasite_info <- trimws(parasite_info)
        
        # Parse parasite names
        if (nchar(parasite_info) > 0) {
          clean_parasite <- gsub("\\?", "", parasite_info)
          parasite_names <- str_extract_all(clean_parasite, "[A-Z][a-z]+\\s+[a-z]+")[[1]]
          
          if (length(parasite_names) > 0) {
            for (parasite_name in parasite_names) {
              parasites_list[[length(parasites_list) + 1]] <- data.frame(
                host = full_species,
                parasite = trimws(parasite_name),
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
      
      # Extract plant interaction information
      pollen_info <- ""
      if (pollen_pos > 0) {
        pollen_start <- pollen_pos + 7  # Skip "Pollen:"
        pollen_end <- nchar(species_text)
        
        # Find where the citation line starts in species_text to exclude it from pollen_info
        if (!is.na(citation_line)) {
          # Try to find the citation by looking for "Genus species Author"
          # This is more reliable than looking for just the genus
          citation_pattern_search <- paste0(current_genus, "\\s+", species_epithet, "\\s+[A-Z][a-z]+")
          citation_in_text <- regexpr(citation_pattern_search, species_text)
          
          # If subspecies, also try with subspecies epithet
          if (citation_in_text < 0 && !is.na(subspecies_epithet)) {
            citation_pattern_search <- paste0(current_genus, "\\s+", subspecies_epithet, "\\s+[A-Z][a-z]+")
            citation_in_text <- regexpr(citation_pattern_search, species_text)
          }
          
          if (citation_in_text > 0 && citation_in_text > pollen_start) {
            pollen_end <- citation_in_text - 1
          }
        }
        
        pollen_info <- substr(species_text, pollen_start, pollen_end)
        pollen_info <- trimws(pollen_info)
      }
      
      # Determine diet breadth
      diet_breadth <- NA
      if (nchar(pollen_info) > 0) {
        
        # 1. Polylectic — check first, to prioritize
        if (grepl("\\b(polylege|polylectic|polytropic)\\b", pollen_info, ignore.case = TRUE)) {
          diet_breadth <- "polylectic"
          
          # 2. Possibly oligolectic — 'possibly' before oligolectic keywords
        } else if (grepl("\\b(possibly|apparently)\\b.{0,20}?\\b(oligolege|oligolectic|oligotrophic)\\b", 
                         pollen_info, ignore.case = TRUE)) {
          diet_breadth <- "possibly oligolectic"
          
          # 3. Oligolectic — plain mention
        } else if (grepl("\\b(oligolege|oligolectic|oligotrophic)\\b", pollen_info, ignore.case = TRUE)) {
          diet_breadth <- "oligolectic"
        }
      }
      
      # Add to species list
      # Validation: Skip if this is a known predator/parasite genus (false positive)
      is_valid_species <- TRUE
      
      # # Common predator/parasite genera that might be mistaken for species
      # predator_genera <- c("Philanthus", "Epeolus", "Triepeolus", "Nomada", "Sphecodes", 
      #                      "Stelis", "Coelioxys", "Melecta")
      # 
      # # Check if the species epithet matches a known predator genus
      # # (This catches cases where "albopilosus" from "Philanthus albopilosus" was detected)
      # if (!is.na(species_epithet)) {
      #   # Get the full species name we created
      #   test_genus <- strsplit(full_species, " ")[[1]][1]
      #   
      #   # If somehow we created a species using a predator genus name, skip it
      #   if (test_genus %in% predator_genera) {
      #     is_valid_species <- FALSE
      #   }
      # }
      
      if (is_valid_species) {
        species_list[[length(species_list) + 1]] <- data.frame(
          higher = full_species,
          author = author_with_year,
          distribution = distribution,
          dietBreadth = diet_breadth,
          speciesGroup = current_species_group,
          stringsAsFactors = FALSE
        )
        cat("  ADDED TO SPECIES LIST (total:", length(species_list), ")\n")
      } else {
        cat("  REJECTED: Failed validation\n")
      }
      
      # IMPROVED INTERACTION CLASSIFICATION
      # Replace the section starting from "# Parse plant interactions" 
      
      # Parse plant interactions
      if (nchar(pollen_info) > 0 && !grepl("^\\s*Unknown\\s*$", pollen_info, ignore.case = TRUE)) {
        
        # --- STEP 1: Split text at visitation boundaries ---
        # Split at phrases that indicate a shift from pollen collection to nectar visitation:
        # - "but visits other flowers"
        # - "but visits these and other flowers"
        # - "but also visits"
        # - "but visits"  
        # - "visiting flowers of"
        
        split_pattern <- "but visits these and other flowers|but visits other flowers|but also visits|but visits|,?\\s*visiting flowers of"
        split_pos <- regexpr(split_pattern, pollen_info, ignore.case = TRUE)[1]
        
        if (split_pos > 0) {
          # Text is split into pollen and nectar sections
          pollen_section <- substr(pollen_info, 1, split_pos - 1)
          nectar_section <- substr(pollen_info, split_pos, nchar(pollen_info))
        } else {
          # No split - determine based on keywords
          if (grepl("collects pollen|obtains pollen|pollen from|oligolege|polylege|especially flowers", pollen_info, ignore.case = TRUE)) {
            pollen_section <- pollen_info
            nectar_section <- ""
          } else {
            pollen_section <- ""
            nectar_section <- pollen_info
          }
        }
        
        # --- STEP 2: Extract plants from each section ---
        extract_plants <- function(text_section) {
          if (nchar(trimws(text_section)) == 0) return(c())
          
          # Clean the text
          clean <- text_section
          clean <- gsub("\\\\\\s*", " ", clean)
          clean <- gsub("\\s+", " ", clean)
          
          # Remove citation patterns
          clean <- gsub(paste0(current_genus, "\\s+[a-z]+\\s+[A-Z][^,]*,\\s*\\d{4}\\..*$"), "", clean)
          clean <- gsub("[A-Z][a-z]+,\\s*\\d{4}\\.\\s+[A-Z].*$", "", clean)
          
          # Remove  narrative phrases that aren't plant data
          clean <- gsub("\\bBased upon[^,.;]+[,.;]?", "", clean, ignore.case = TRUE)
          clean <- gsub("\\bthis species evidently\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("\\bthe mouth parts of the female,?\\s+", "", clean, ignore.case = TRUE)
          
          # Remove descriptive prefixes
          clean <- gsub("Unknown,?\\s+but\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("^\\s*(Apparently|Presumably|Possibly|Probably)\\s+(an\\s+)?oligolege\\s+(of|from)?\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("^\\s*(Apparently|Presumably|Possibly|Probably)\\s+(an\\s+)?polylege\\s+(of|from)?\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("^\\s*(Apparently|Presumably|Possibly|Probably|Evidently)\\s+(collects pollen|obtains pollen)\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("^\\s*(collects pollen|obtains pollen|pollen)\\s+(from|of)\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("^\\s*(Oligolege|Polylege|Polylectic)\\s+(of|from)?\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("^\\s*(visits|visiting)\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("^\\s*(but visits|but also visits)\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub(".*?flowers\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("^(of|including|presumably)\\s+", "", clean, ignore.case = TRUE)
          
          # Handle "including" and common typos
          clean <- gsub("\\s+ineluding\\s+", ", ", clean, ignore.case = TRUE)  # typo in source
          clean <- gsub("\\s+(including|such as|e\\.g\\.,?)\\s+(but not restricted to\\s+)?", ", ", clean, ignore.case = TRUE)
          clean <- gsub("^\\s*(autumnal|vernal|spring|fall)\\s+flowering\\s+", "", clean, ignore.case = TRUE)
          clean <- gsub("\\s*for nectar\\b", "", clean, ignore.case = TRUE)
          clean <- gsub("\\s+presumably\\b", "", clean, ignore.case = TRUE)
          clean <- gsub("\\s+primarily from\\s+", ", ", clean, ignore.case = TRUE)
          clean <- gsub(",?\\s+visiting flowers of\\s+", ", ", clean, ignore.case = TRUE)
          
          clean <- gsub("\\s+the flowers of\\s+", " ", clean, ignore.case = TRUE)
          
          clean <- trimws(clean)
          
          # Split by delimiters
          parts <- unlist(strsplit(clean, "[,;]|\\sand\\s|\\sor\\s"))
          
          plant_list <- c()
          last_genus <- NULL
          
          for (part in parts) {
            part <- trimws(part)
            if (nchar(part) < 2) next
            
            # Full binomial (Genus species)
            binomial <- str_match(part, "^([A-Z][a-z]+)\\s+([a-z]+)")
            if (!is.na(binomial[1])) {
              plant_name <- paste(binomial[2], binomial[3])
              plant_list <- c(plant_list, plant_name)
              last_genus <- binomial[2]
              next
            }
            
            # Abbreviated with period (G. species)
            abbrev <- str_match(part, "^([A-Z])\\.\\s+([a-z]+)")
            if (!is.na(abbrev[1]) && !is.null(last_genus)) {
              plant_name <- paste(last_genus, abbrev[3])
              plant_list <- c(plant_list, plant_name)
              next
            }
            
            # Abbreviated without period (G species)
            abbrev_no_period <- str_match(part, "^([A-Z])\\s+([a-z]+)")
            if (!is.na(abbrev_no_period[1]) && !is.null(last_genus)) {
              plant_name <- paste(last_genus, abbrev_no_period[3])
              plant_list <- c(plant_list, plant_name)
              next
            }
            
            # Family name (-aceae or -ae ending)
            family_match <- str_match(part, "^([A-Z][a-z]+a[ce]ae?)\\b")
            if (!is.na(family_match[1])) {
              plant_list <- c(plant_list, family_match[2])
              # Don't update last_genus for families
              next
            }
            
            # Just genus
            genus_only <- str_match(part, "([A-Z][a-z]+)")
            if (!is.na(genus_only[1])) {
              genus_name <- genus_only[2]
              
              # Exclude non-plant terms
              exclude_terms <- c("Unknown", "Oligolege", "Polylege", "Presumably", "Parasite", "Predator", 
                                 "Biology", "Ecology", "Monotypic", "Mexico", "California", "Arizona",
                                 "Polylectic", "Polytropic", "Males", "Females", "Spring",
                                 "Fall", "Autumnal", "Vernal", "Pollen", "Texas", "Colorado",
                                 "Apparently", "Possibly", "Various", "Typically", "Generally",
                                 "Nevada", "Utah", "Idaho", "Wyoming", "Montana", "Alberta",
                                 "British", "Columbia", "New", "North", "South", "East", "West",
                                 "Produces", "The", "Of", "In", "At", "To", "From", "For",
                                 "With", "And", "Or", "But", "Also", "Only", "Primarily",
                                 "Ann", "Mag", "Proc", "Trans", "Bull", "Bul", "Jour", "Rev",
                                 "Ent", "Soc", "Nat", "Hist", "Sci", "Acad", "Amer", "Zool",
                                 "Univ", "Dept", "Tech", "Agr", "Expt", "Sta", "Fig", "Figs",
                                 "Vol", "Contrib", "Smithsn", "Can", "Canad", "Calif",
                                 "Including", "Visiting", "Especially", "Based", "Evidently",
                                 "Probably", "Collects", "Upon", "Mouth", "Parts", "Female",
                                 "Species", "Ineluding")  # last one is typo in source
              
              if (!genus_name %in% exclude_terms && nchar(genus_name) > 2) {
                plant_list <- c(plant_list, genus_name)
                last_genus <- genus_name
              }
            }
          }
          
          return(unique(plant_list))
        }
        
        # Extract plants from each section
        pollen_plants <- extract_plants(pollen_section)
        nectar_plants <- extract_plants(nectar_section)
        
        # --- STEP 3: Create interaction records ---
        # Pollen plants get collectsPollenOf
        for (plant in pollen_plants) {
          interactions_list[[length(interactions_list) + 1]] <- data.frame(
            higher = full_species,
            interactionTypeName = "collectsPollenOf",
            lower = trimws(plant),
            stringsAsFactors = FALSE
          )
        }
        
        # Nectar plants get visitsFlowersOf
        for (plant in nectar_plants) {
          # Don't add duplicates if plant appears in both sections
          if (!(plant %in% pollen_plants)) {
            interactions_list[[length(interactions_list) + 1]] <- data.frame(
              higher = full_species,
              interactionTypeName = "visitsFlowersOf",
              lower = trimws(plant),
              stringsAsFactors = FALSE
            )
          }
        }
      }
      next
    }
    
    i <- i + 1
  }
  
  # Convert lists to data frames
  interactions_df <- bind_rows(interactions_list)
  species_df <- bind_rows(species_list)
  parasites_df <- bind_rows(parasites_list)
  
  # Remove duplicates
  if (nrow(interactions_df) > 0) {
    interactions_df <- interactions_df %>% distinct()
  }
  if (nrow(species_df) > 0) {
    species_df <- species_df %>% distinct()
  }
  if (nrow(parasites_df) > 0) {
    parasites_df <- parasites_df %>% distinct()
  }
  
  return(list(
    interactions = interactions_df,
    species = species_df,
    parasites = parasites_df
  ))
}

#### Parse file ####
results <- parse_hymenoptera_book("Data_Raw/Interactions/BigBook/BigBookOfHymenoptera_Apoidea.txt")

interactions <- results$interactions
# species <- results$species
## parasites list unreliable right now
# parasites <- results$parasites


# fix capitalization of genus names
fix_genus_case <- function(x) {
  sapply(strsplit(x, "\\s+"), function(parts) {
    if (length(parts) == 0) return(NA_character_)
    parts[1] <- str_to_title(parts[1])
    paste(parts, collapse = " ")
  })
}

interactions <- interactions %>%
  mutate(higher = fix_genus_case(higher))

#### Save raw interactions data scraped ####
write.csv(interactions, "Data_Raw/Interactions/BigBook/bigBook_scraped_interactions.csv")