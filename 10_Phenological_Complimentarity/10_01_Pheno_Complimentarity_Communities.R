#### Load Packages ####
library(dplyr)
library(tidyr)
library(parallel)
library(foreach)
library(doParallel)

#### Load Data ####
plantSpp = read.csv("Data_Calscape/Interactions/calscape_interactions_allPredicted.csv") %>%
  pull(lower) %>%
  unique()

## Load phenometrics for calscape plants
# load phenology data
phenometrics <- read.csv("Data_Clean/Phenology/phenology_best_available.csv") %>%
  # filter to only plants 
  filter(genus_species %in% plantSpp) %>%
  dplyr::select(genus_species, JEPCODE, phen_start, phen_end) %>%
  rename(scientificName = genus_species,
         phen_05 = phen_start,
         phen_95 = phen_end) %>%
  filter(!is.na(phen_05) & !is.na(phen_95)) %>%
  unique()


# Calflora plant community information
calflora_communities <- read.csv("Data_Clean/Species_Traits_Attributes/Calflora_Communities.csv") %>%
  dplyr::select(genus_species_harm, community) %>%
  rename(lower = genus_species_harm) %>%
  unique()

calflora_communities = na.omit(calflora_communities)

## Make ecoregions list
ecoregions = unique(phenometrics$JEPCODE)

# Set up the parallel backend
cl <- makeCluster(24)
registerDoParallel(cl)

# Function to find plants that share at least one community with the focal plant
find_community_sharing_plants <- function(focal_plant, community_matrix, ecoregion_plants) {
  # Check if focal plant is in community data
  if (!(focal_plant %in% rownames(community_matrix))) {
    # If plant has no community data, return all plants in the ecoregion except itself
    return(setdiff(ecoregion_plants, focal_plant))
  }
  
  # if it does have community information, filter for only other plants in those communities
  focal_communities <- community_matrix[focal_plant, ]
  # Find plants that share at least one community
  shared_community_plants <- character(0)
  
  for (plant in rownames(community_matrix)) {
    if (plant == focal_plant) next
    
    plant_communities <- community_matrix[plant, ]
    # Check if there's an overlap (at least one shared community)
    if (any(focal_communities & plant_communities)) {
      shared_community_plants <- c(shared_community_plants, plant)
    }
  }
  
  return(shared_community_plants)
}

# Function to find 10 best companion plants for each plant in an ecoregion,
# constrained to plants that share communities
find_community_constrained_complimentary_phenology <- function(ecoregion) {
  # Filter data for the specific ecoregion
  ecoregion_data <- phenometrics %>%
    filter(JEPCODE %in% c(ecoregion))
  
  # Get unique plants in this ecoregion
  ecoregion_plants <- unique(ecoregion_data$scientificName)
  
  # Create a daily binary phenomatrix using the phenmetrics in ecoregion_data
  # Create an empty phenology matrix first
  pheno_matrix <- matrix(F, 
                         nrow = length(ecoregion_plants), 
                         ncol = 365,
                         dimnames = list(ecoregion_plants, 1:365))
  
  # Now fill in the flowering days with 1s
  for (i in 1:nrow(ecoregion_data)) {
    species <- ecoregion_data$scientificName[i]
    start_day <- max(1, min(365, ceiling(ecoregion_data$phen_05[i])))
    end_day   <- max(1, min(365, floor(ecoregion_data$phen_95[i])))
    
    # Set flowering days to 1 using indices
    species_idx <- which(rownames(pheno_matrix) == species)
    pheno_matrix[species_idx, start_day:end_day] <- T
  }
  
  # Filter community data to include only plants in this ecoregion
  filtered_communities <- calflora_communities %>%
    filter(lower %in% ecoregion_plants)
  
  # Generate a binary plant x community matrix
  plant_community_matrix <- filtered_communities %>%
    # Assuming the community column name is 'community'
    # Adjust column name if different
    dplyr::select(lower, community) %>%
    unique() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = community,
                values_from = present,
                values_fill = 0)
  
  # Convert to matrix format for efficient operations
  plant_names <- plant_community_matrix$lower
  community_matrix <- as.matrix(plant_community_matrix[, -1])
  rownames(community_matrix) <- plant_names
  
  # Initialize result dataframe
  companion_results <- data.frame(
    focal_plant = character(0),
    companion_plant = character(0),
    new_days = integer(0),
    total_days = integer(0),
    rank = integer(0),
    ecoregion = character(0)
  )
  
  # For each plant in the ecoregion, find companions
  for (focal_plant in ecoregion_plants) {
    # Skip if focal plant isn't in pheno_matrix (shouldn't happen but just in case)
    if (!(focal_plant %in% rownames(pheno_matrix))) next
    
    # Get pheno vector for the focal plant
    focal_pheno_vector <- pheno_matrix[focal_plant, ]
    
    # Get plants that share communities with the focal plant
    community_sharing_plants <- find_community_sharing_plants(focal_plant, community_matrix, ecoregion_plants)
    
    # If no plants share communities (shouldn't happen with the updated function, but just in case)
    if (length(community_sharing_plants) == 0) {
      # Fall back to using all plants in the ecoregion except the focal plant
      community_sharing_plants <- setdiff(ecoregion_plants, focal_plant)
    }
    
    # Filter to only include plants that are in the pheno_matrix
    candidate_plants <- intersect(community_sharing_plants, rownames(pheno_matrix))
    
    # If no candidates after filtering, skip
    if (length(candidate_plants) == 0) next
    
    # Calculate the complementarity of each candidate plant with the focal plant only
    # (1-1 comparison rather than iterative selection)
    new_days_counts <- sapply(candidate_plants, function(plant) {
      plant_days <- pheno_matrix[plant, ]
      # Calculate how many new days this plant adds to the focal plant's phenology
      sum(plant_days & !focal_pheno_vector)
    })
    
    # Calculate total days covered by combining focal plant with each candidate
    total_days_counts <- sapply(candidate_plants, function(plant) {
      plant_days <- pheno_matrix[plant, ]
      # Sum of days covered by either focal plant or the candidate plant
      sum(focal_pheno_vector | plant_days)
    })
    
    # Create a data frame with all candidates and their metrics
    candidates_df <- data.frame(
      plant = candidate_plants,
      new_days = new_days_counts,
      total_days = total_days_counts
    )
    
    # Sort by new_days in descending order
    candidates_df <- candidates_df[order(-candidates_df$new_days), ]
    
    # Filter to remove plants that dont add new days
    candidates_df <- candidates_df[candidates_df$new_days > 0, ]
    
    # Select the top 10 plants (or fewer if there aren't 10 candidates)
    top_plants <- head(candidates_df, 20)
    
    # Add to results dataframe
    if (nrow(top_plants) > 0) {
      companion_results <- rbind(companion_results, data.frame(
        focal_plant     = focal_plant,
        companion_plant = top_plants$plant,
        new_days        = top_plants$new_days,
        total_days      = top_plants$total_days,
        rank            = seq_len(nrow(top_plants)),
        JEPCODE         = ecoregion
      ))
    }
  }
  
  return(companion_results)
}

## Process all ecoregions in parallel
complimentarity_results <- foreach(ecoregion = ecoregions, .packages = c("dplyr", "tidyr"), .combine = rbind) %dopar% {
  find_community_constrained_complimentary_phenology(ecoregion)
}


#### State level pheno complimentarity ####
# Same algorithm but pooled across all ecoregions, with candidates filtered to plants co-occurring in >= 50% of the focal plant's ecoregions.

## Load and combine state-level phenology (Type 1 + Type 2 + Type 3)
phen_simple <- read.csv("Data_Clean/Phenology/phenology_simple.csv")
phen_type2  <- read.csv("Data_Clean/Phenology/phenology_type2.csv")
phen_type3  <- read.csv("Data_Clean/Phenology/phenology_type3.csv")

# Prepare
state_phenology <- bind_rows(phen_simple, phen_type2, phen_type3) %>%
  filter(genus_species %in% plantSpp) %>%
  dplyr::select(genus_species, phen_start, phen_end) %>%
  rename(scientificName = genus_species, phen_05 = phen_start, phen_95 = phen_end) %>%
  filter(!is.na(phen_05) & !is.na(phen_95)) %>%
  unique()

## Build state-level pheno_matrix
state_plants <- unique(state_phenology$scientificName)

state_pheno_matrix <- matrix(FALSE,
                             nrow = length(state_plants),
                             ncol = 365,
                             dimnames = list(state_plants, 1:365))

for (i in 1:nrow(state_phenology)) {
  species <- state_phenology$scientificName[i]
  start_day <- max(1, min(365, ceiling(state_phenology$phen_05[i])))
  end_day   <- max(1, min(365, floor(state_phenology$phen_95[i])))
  species_idx <- which(rownames(state_pheno_matrix) == species)
  state_pheno_matrix[species_idx, start_day:end_day] <- TRUE
}

## Build state-level community matrix
state_filtered_communities <- calflora_communities %>%
  filter(lower %in% state_plants)

state_plant_community_df <- state_filtered_communities %>%
  dplyr::select(lower, community) %>%
  unique() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = community,
              values_from = present,
              values_fill = 0)

state_comm_names <- state_plant_community_df$lower
state_community_matrix <- as.matrix(state_plant_community_df[, -1])
rownames(state_community_matrix) <- state_comm_names

## Build plant x ecoregion co-occurrence matrix (from ecoregion-level phenometrics)
plant_eco_df <- phenometrics %>%
  dplyr::select(scientificName, JEPCODE) %>%
  unique() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = JEPCODE,
              values_from = present,
              values_fill = 0)

eco_plant_names <- plant_eco_df$scientificName
state_eco_matrix <- as.matrix(plant_eco_df[, -1])
rownames(state_eco_matrix) <- eco_plant_names

## Processing a single focal plant at the state level
process_state_focal_plant_phenology <- function(focal_plant,
                                                state_pheno_matrix,
                                                state_community_matrix,
                                                state_eco_matrix,
                                                all_plants) {
  
  # Empty result template
  empty_result <- data.frame(focal_plant = character(0), companion_plant = character(0),
                             new_days = integer(0), total_days = integer(0),
                             rank = integer(0), JEPCODE = character(0))
  
  if (!(focal_plant %in% rownames(state_pheno_matrix))) return(empty_result)
  
  focal_pheno_vector <- state_pheno_matrix[focal_plant, ]
  
  ## Community filter (identical to ecoregion level)
  community_sharing <- find_community_sharing_plants(
    focal_plant, state_community_matrix, all_plants
  )
  if (length(community_sharing) == 0) {
    community_sharing <- setdiff(all_plants, focal_plant)
  }
  
  ## Co-occurrence filter: companion must share >= 50% of focal's ecoregions
  if (focal_plant %in% rownames(state_eco_matrix)) {
    focal_eco     <- state_eco_matrix[focal_plant, ]
    n_focal_eco   <- sum(focal_eco)
    shared_counts <- as.vector(state_eco_matrix %*% focal_eco)
    names(shared_counts) <- rownames(state_eco_matrix)
    cooccurring <- names(shared_counts)[shared_counts >= n_focal_eco * 0.5 &
                                          names(shared_counts) != focal_plant]
  } else {
    cooccurring <- setdiff(all_plants, focal_plant)
  }
  
  ## Intersect both filters, then restrict to pheno_matrix rows
  candidate_plants <- intersect(community_sharing, cooccurring)
  candidate_plants <- intersect(candidate_plants, rownames(state_pheno_matrix))
  
  if (length(candidate_plants) == 0) return(empty_result)
  
  ## Non-greedy complementarity (same algorithm as ecoregion level)
  # Calculate new days each candidate adds relative to focal plant only
  new_days_counts <- sapply(candidate_plants, function(plant) {
    sum(state_pheno_matrix[plant, ] & !focal_pheno_vector)
  })
  
  total_days_counts <- sapply(candidate_plants, function(plant) {
    sum(focal_pheno_vector | state_pheno_matrix[plant, ])
  })
  
  candidates_df <- data.frame(
    plant = candidate_plants,
    new_days = new_days_counts,
    total_days = total_days_counts
  )
  
  # which plants add the most new days pheno (break ties by sorting by all_days)
  candidates_df <- candidates_df[order(-candidates_df$new_days, -candidates_df$total_days), ]
  
  # remove days that don't add anything new
  candidates_df <- candidates_df[candidates_df$new_days > 0, ]
  
  top_plants <- head(candidates_df, 20)
  
  results <- data.frame(
    focal_plant     = rep(focal_plant, nrow(top_plants)),
    companion_plant = top_plants$plant,
    new_days        = top_plants$new_days,
    total_days      = top_plants$total_days,
    rank            = 1:nrow(top_plants),
    JEPCODE         = "California"
  )
  
  return(results)
}

## Run state-level phenological complementarity in parallel
state_complimentarity_results <- foreach(
  fp = state_plants,
  .packages = c("dplyr", "tidyr"),
  .combine  = rbind,
  .export   = c("find_community_sharing_plants", "process_state_focal_plant_phenology",
                "state_pheno_matrix", "state_community_matrix", "state_eco_matrix",
                "state_plants")
) %dopar% {
  process_state_focal_plant_phenology(
    focal_plant            = fp,
    state_pheno_matrix     = state_pheno_matrix,
    state_community_matrix = state_community_matrix,
    state_eco_matrix       = state_eco_matrix,
    all_plants             = state_plants
  )
}

#### Combine ecoregion-level and state-level results and save ####
complimentarity_results <- rbind(complimentarity_results, state_complimentarity_results)

# Stop the cluster
stopCluster(cl)

# Write detailed results to CSV
write.csv(complimentarity_results, "Data_Calscape/Companion_Plants/community_constrained_phenological_complimentarity.csv", row.names = FALSE)

