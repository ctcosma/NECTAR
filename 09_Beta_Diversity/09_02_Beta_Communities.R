#### Load packages ####
library(dplyr)
library(tidyr)
library(parallel)
library(foreach)
library(doParallel)

#### Load and clean data ####
## Load interaction data 
interaction_data = read.csv("Data_Calscape/Interactions/calscape_interactions_allPredicted.csv")

## Calflora plant community information
calflora_communities <- read.csv("Data_Clean/Species_Traits_Attributes/Calflora_Communities.csv") %>%
  dplyr::select(genus_species_harm, community) %>%
  rename(lower = genus_species_harm) %>%
  unique()

calflora_communities = na.omit(calflora_communities)

## Interaction cleaning
interaction_data_nectar <- interaction_data %>%
  # Beta diversity for only nectar interactions
  filter(interactionTypeName %in% c("visitsFlowersOf", "collectsPollenOf_general", "collectsPollenOf_specialist")) %>%
  # select clean columns
  dplyr::select(higher, lower, JEPCODE) %>%
  mutate(potentialOrConfirmed = 1) %>%
  unique()

interaction_data_host <- interaction_data %>%
  # Beta diversity for only host interactions
  filter(interactionTypeName == "hasHost") %>%
  # select clean columns
  dplyr::select(higher, lower, JEPCODE) %>%
  mutate(potentialOrConfirmed = 1) %>%
  unique()

## Make ecoregions list
ecoregions = unique(interaction_data$JEPCODE)

# Set up the parallel backend
cl <- makeCluster(24)
registerDoParallel(cl)

#### Ecoregion companions ####
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

# Function to find 20 best companion plants for each plant in an ecoregion,
# constrained to plants that share communities
find_community_constrained_companions <- function(ecoregion, interaction_data_subset) {
  # Filter data for the specific ecoregion
  ecoregion_data <- interaction_data_subset %>%
    filter(JEPCODE %in% c(ecoregion))
  
  # Get unique plants and pollinators in this ecoregion
  ecoregion_plants <- unique(ecoregion_data$lower)
  pollinators <- unique(ecoregion_data$higher)
  
  # Create a plant-pollinator interaction matrix
  interaction_matrix <- ecoregion_data %>%
    select(lower, higher, potentialOrConfirmed) %>%
    pivot_wider(names_from = higher, 
                values_from = potentialOrConfirmed,
                values_fill = 0) %>%
    tibble::column_to_rownames("lower")
  
  # Convert to logical matrix for more efficient operations
  interaction_matrix <- as.matrix(interaction_matrix) > 0
  
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
    added_pollinators = integer(0),
    total_pollinators = integer(0),
    rank = integer(0),
    ecoregion = character(0)
  )
  
  # For each plant in the ecoregion, find companions
  for (focal_plant in ecoregion_plants) {
    # Skip if focal plant isn't in interaction_matrix (shouldn't happen but just in case)
    if (!(focal_plant %in% rownames(interaction_matrix))) next
    
    # Get pollinators of the focal plant
    focal_pollinators <- interaction_matrix[focal_plant, ]
    
    # Get plants that share communities with the focal plant
    community_sharing_plants <- find_community_sharing_plants(focal_plant, community_matrix, ecoregion_plants)
    
    # If no plants share communities (shouldn't happen with the updated function, but just in case)
    if (length(community_sharing_plants) == 0) {
      # Fall back to using all plants in the ecoregion except the focal plant
      community_sharing_plants <- setdiff(ecoregion_plants, focal_plant)
    }
    
    # Filter to only include plants that are in the interaction_matrix
    candidate_plants <- intersect(community_sharing_plants, rownames(interaction_matrix))
    
    # If no candidates after filtering, skip
    if (length(candidate_plants) == 0) next
    
    # Initialize for this focal plant
    covered_pollinators <- focal_pollinators
    selected_companions <- character(0)
    
    # Select up to 10 companion plants
    for (rank in 1:20) {
      if (length(candidate_plants) == 0) break
      
      # Calculate the number of new pollinators each plant would add
      new_pollinator_counts <- sapply(candidate_plants, function(plant) {
        plant_pollinators <- interaction_matrix[plant, ]
        sum(plant_pollinators & !covered_pollinators)
      })
      
      # If no new pollinators can be added, stop
      if (max(new_pollinator_counts) == 0) break
      
      # Find plant with max new pollinators
      ### TEAGAN CHECK: Check how ties are being handled. Maybe should break ties in new pollinators by total pollinators supported by each plant?
      best_plant_index <- which.max(new_pollinator_counts)
      best_plant <- candidate_plants[best_plant_index]
      best_count <- new_pollinator_counts[best_plant_index]
      
      # Update covered pollinators
      plant_pollinators <- interaction_matrix[best_plant, ]
      covered_pollinators <- covered_pollinators | plant_pollinators
      
      # Add best plant to selected list
      selected_companions <- c(selected_companions, best_plant)
      
      # Remove selected plant from candidate list
      candidate_plants <- setdiff(candidate_plants, best_plant)
      
      # Record this companion plant
      companion_results <- rbind(companion_results, data.frame(
        focal_plant = focal_plant,
        companion_plant = best_plant,
        added_pollinators = best_count,
        total_pollinators = sum(covered_pollinators),
        rank = rank,
        JEPCODE = ecoregion
      ))
    }
  }
  
  return(companion_results)
}

## Process all ecoregions in parallel
# nectar plants
all_companion_nectar_results <- foreach(ecoregion = ecoregions, .packages = c("dplyr", "tidyr"), .combine = rbind) %dopar% {
  find_community_constrained_companions(ecoregion, interaction_data_nectar)
}
# host plants
all_companion_host_results <- foreach(ecoregion = ecoregions, .packages = c("dplyr", "tidyr"), .combine = rbind) %dopar% {
  find_community_constrained_companions(ecoregion, interaction_data_host)
}

#### state level companion plants #####
# Same greedy algorithm but pooled across all ecoregions, with candidates filtered to plants co-occurring in >= 50% of the focal plant's ecoregions.

# prepare state-level data structures from an interaction dataset
prepare_state_data <- function(interaction_data_subset) {
  
  # 1. State-level interaction matrix (unique lower–higher pairs, no ecoregion)
  state_data <- interaction_data_subset %>%
    dplyr::select(higher, lower, potentialOrConfirmed) %>%
    unique()
  
  all_plants <- unique(state_data$lower)
  
  state_interaction_matrix <- state_data %>%
    dplyr::select(lower, higher, potentialOrConfirmed) %>%
    pivot_wider(names_from = higher,
                values_from = potentialOrConfirmed,
                values_fill = 0) %>%
    tibble::column_to_rownames("lower")
  state_interaction_matrix <- as.matrix(state_interaction_matrix) > 0
  
  # 2. State-level community matrix (same logic as ecoregion level)
  filtered_communities <- calflora_communities %>%
    filter(lower %in% all_plants)
  
  plant_community_df <- filtered_communities %>%
    dplyr::select(lower, community) %>%
    unique() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = community,
                values_from = present,
                values_fill = 0)
  
  comm_plant_names <- plant_community_df$lower
  state_community_matrix <- as.matrix(plant_community_df[, -1])
  rownames(state_community_matrix) <- comm_plant_names
  
  # 3. Plant × ecoregion binary matrix (for co-occurrence filtering)
  plant_eco_df <- interaction_data_subset %>%
    dplyr::select(lower, JEPCODE) %>%
    unique() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = JEPCODE,
                values_from = present,
                values_fill = 0)
  
  eco_plant_names <- plant_eco_df$lower
  state_eco_matrix <- as.matrix(plant_eco_df[, -1])
  rownames(state_eco_matrix) <- eco_plant_names
  
  return(list(
    interaction_matrix = state_interaction_matrix,
    community_matrix  = state_community_matrix,
    eco_matrix        = state_eco_matrix,
    all_plants        = all_plants
  ))
}

# get companions for single plant state level
process_state_focal_plant <- function(focal_plant,
                                      state_interaction_matrix,
                                      state_community_matrix,
                                      state_eco_matrix,
                                      all_plants) {
  
  if (!(focal_plant %in% rownames(state_interaction_matrix))) {
    return(data.frame(focal_plant = character(0), companion_plant = character(0),
                      added_pollinators = integer(0), total_pollinators = integer(0),
                      rank = integer(0), JEPCODE = character(0)))
  }
  
  focal_pollinators <- state_interaction_matrix[focal_plant, ]
  
  # Community filter (identical to ecoregion level)
  community_sharing <- find_community_sharing_plants(
    focal_plant, state_community_matrix, all_plants
  )
  if (length(community_sharing) == 0) {
    community_sharing <- setdiff(all_plants, focal_plant)
  }
  
  # Co-occurrence filter: companion must share >= 50% of focal's ecoregions
  if (focal_plant %in% rownames(state_eco_matrix)) {
    focal_eco       <- state_eco_matrix[focal_plant, ]
    n_focal_eco     <- sum(focal_eco)
    shared_counts   <- as.vector(state_eco_matrix %*% focal_eco)
    names(shared_counts) <- rownames(state_eco_matrix)
    cooccurring <- names(shared_counts)[shared_counts >= n_focal_eco * 0.5 &
                                          names(shared_counts) != focal_plant]
  } else {
    # No ecoregion info — fall back to all plants
    cooccurring <- setdiff(all_plants, focal_plant)
  }
  
  # Intersect both filters, then restrict to interaction matrix rows
  candidate_plants <- intersect(community_sharing, cooccurring)
  candidate_plants <- intersect(candidate_plants, rownames(state_interaction_matrix))
  
  if (length(candidate_plants) == 0) {
    return(data.frame(focal_plant = character(0), companion_plant = character(0),
                      added_pollinators = integer(0), total_pollinators = integer(0),
                      rank = integer(0), JEPCODE = character(0)))
  }
  
  # Greedy companion selection
  results <- data.frame(focal_plant = character(0), companion_plant = character(0),
                        added_pollinators = integer(0), total_pollinators = integer(0),
                        rank = integer(0), JEPCODE = character(0))
  
  covered_pollinators <- focal_pollinators
  
  for (r in 1:20) {
    if (length(candidate_plants) == 0) break
    
    new_pollinator_counts <- sapply(candidate_plants, function(plant) {
      sum(state_interaction_matrix[plant, ] & !covered_pollinators)
    })
    
    if (max(new_pollinator_counts) == 0) break
    
    best_idx   <- which.max(new_pollinator_counts)
    best_plant <- candidate_plants[best_idx]
    best_count <- new_pollinator_counts[best_idx]
    
    covered_pollinators <- covered_pollinators | state_interaction_matrix[best_plant, ]
    candidate_plants    <- setdiff(candidate_plants, best_plant)
    
    results <- rbind(results, data.frame(
      focal_plant     = focal_plant,
      companion_plant = best_plant,
      added_pollinators = best_count,
      total_pollinators = sum(covered_pollinators),
      rank   = r,
      JEPCODE = "California"
    ))
  }
  
  return(results)
}

#### Run state-level companions for NECTAR interactions ####
state_nectar_data <- prepare_state_data(interaction_data_nectar)

state_companion_nectar_results <- foreach(
  fp = state_nectar_data$all_plants,
  .packages = c("dplyr", "tidyr"),
  .combine  = rbind,
  .export   = c("find_community_sharing_plants", "process_state_focal_plant",
                "state_nectar_data")
) %dopar% {
  process_state_focal_plant(
    focal_plant             = fp,
    state_interaction_matrix = state_nectar_data$interaction_matrix,
    state_community_matrix  = state_nectar_data$community_matrix,
    state_eco_matrix        = state_nectar_data$eco_matrix,
    all_plants              = state_nectar_data$all_plants
  )
}

#### Run state-level companions for HOST interactions ####
state_host_data <- prepare_state_data(interaction_data_host)

state_companion_host_results <- foreach(
  fp = state_host_data$all_plants,
  .packages = c("dplyr", "tidyr"),
  .combine  = rbind,
  .export   = c("find_community_sharing_plants", "process_state_focal_plant",
                "state_host_data")
) %dopar% {
  process_state_focal_plant(
    focal_plant             = fp,
    state_interaction_matrix = state_host_data$interaction_matrix,
    state_community_matrix  = state_host_data$community_matrix,
    state_eco_matrix        = state_host_data$eco_matrix,
    all_plants              = state_host_data$all_plants
  )
}

#### Combine ecoregion and state results and save ####
all_companion_nectar_results <- rbind(all_companion_nectar_results, state_companion_nectar_results)
all_companion_host_results   <- rbind(all_companion_host_results, state_companion_host_results)

# Stop the cluster
stopCluster(cl)

# Write detailed results to CSV
write.csv(all_companion_nectar_results, "Data_Calscape/Companion_Plants/community_constrained_companions_nectar.csv", row.names = FALSE)
write.csv(all_companion_host_results, "Data_Calscape/Companion_Plants/community_constrained_companions_host.csv", row.names = FALSE)