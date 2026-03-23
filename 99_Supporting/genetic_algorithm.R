#### Genetic Algorithm Helper Functions ####
# Dependencies: dplyr, tidyr, tibble
# Functions:
# - run_genetic_algorithm:     Core GA engine (takes a matrix, returns results)
# - build_interaction_matrix:  Builds a logical plant x pollinator matrix
# - find_complementary_plants: Wrapper to find complementary plant set, simple outputs
# - find_top_mix:              Wrapper to find complementary plant set (with taxon summaries,
#                              filtered plant/pollinator subsets, fitness history)

#### Core GA ####
# Finds the plant combination that maximizes unique pollinators covered.
# Returns: best_solution, best_fitness, best_solutions (all tied), fitness_history.
run_genetic_algorithm <- function(interaction_matrix,
                                  n_plants,
                                  n_generations = 100,
                                  population_size = 200,
                                  mutation_rate = 0.1) {
  
  all_plants <- rownames(interaction_matrix)
  
  # Fitness function: count unique pollinators covered by a plant combo
  calculate_fitness <- function(plant_combo) {
    covered <- apply(interaction_matrix[plant_combo, , drop = FALSE], 2, any)
    sum(covered)
  }
  
  # Unique key for a plant combination (for deduplication of tied solutions)
  combo_key <- function(plant_combo) {
    paste(sort(plant_combo), collapse = "|||")
  }
  
  # Initialize population with random plant combinations
  population <- replicate(population_size,
                          sample(all_plants, n_plants),
                          simplify = FALSE)
  
  best_fitness <- -Inf
  best_solutions <- list()
  fitness_history <- numeric(n_generations)
  
  for (gen in 1:n_generations) {
    fitness_scores <- sapply(population, calculate_fitness)
    
    current_max <- max(fitness_scores)
    
    # Track all solutions tied at the best fitness
    if (current_max > best_fitness) {
      best_fitness <- current_max
      best_solutions <- list()
      max_indices <- which(fitness_scores == current_max)
      for (idx in max_indices) {
        key <- combo_key(population[[idx]])
        best_solutions[[key]] <- population[[idx]]
      }
    } else if (current_max == best_fitness) {
      max_indices <- which(fitness_scores == current_max)
      for (idx in max_indices) {
        key <- combo_key(population[[idx]])
        if (!key %in% names(best_solutions)) {
          best_solutions[[key]] <- population[[idx]]
        }
      }
    }
    
    fitness_history[gen] <- best_fitness
    
    # Tournament selection
    select_parent <- function() {
      tournament_size <- 5
      tournament_indices <- sample(length(population), tournament_size)
      tournament_fitness <- fitness_scores[tournament_indices]
      population[[tournament_indices[which.max(tournament_fitness)]]]
    }
    
    # Next generation
    new_population <- list()
    
    # Elitism: keep best 10%
    n_elite <- max(1, floor(population_size * 0.1))
    elite_indices <- order(fitness_scores, decreasing = TRUE)[1:n_elite]
    new_population[1:n_elite] <- population[elite_indices]
    
    # Fill rest through crossover and mutation
    for (i in (n_elite + 1):population_size) {
      parent1 <- select_parent()
      parent2 <- select_parent()
      
      # Crossover
      n_from_parent1 <- sample(1:(n_plants - 1), 1)
      child <- c(sample(parent1, n_from_parent1),
                 sample(parent2, n_plants - n_from_parent1))
      
      # Ensure no duplicates
      child <- unique(child)
      while (length(child) < n_plants) {
        child <- c(child, sample(setdiff(all_plants, child), 1))
      }
      child <- child[1:n_plants]
      
      # Mutation
      if (runif(1) < mutation_rate) {
        n_mutations <- sample(1:2, 1)
        positions <- sample(n_plants, n_mutations)
        child[positions] <- sample(setdiff(all_plants, child), n_mutations)
      }
      
      new_population[[i]] <- child
    }
    
    population <- new_population
    
    if (gen %% 10 == 0) {
      message(paste("Generation", gen, "- Best fitness:", best_fitness,
                    "- Unique solutions found:", length(best_solutions)))
    }
  }
  
  message(paste("Best combination(s) support", best_fitness, "pollinator species"))
  message(paste("Found", length(best_solutions), "unique combination(s) with this score"))
  
  list(
    best_solution  = best_solutions[[1]],
    best_fitness   = best_fitness,
    best_solutions = best_solutions,
    fitness_history = fitness_history
  )
}

#### Logical interaction matrix ####
# Build a logical plant x pollinator interaction matrix from a dataframe,
# filtered to the provided plant and pollinator vectors.
build_interaction_matrix <- function(interactions_df,
                                     plants,
                                     pollinators,
                                     plant_col = "targetTaxonName_harm",
                                     pollinator_col = "sourceTaxonName_harm") {
  
  mat <- interactions_df %>%
    dplyr::select(all_of(c(plant_col, pollinator_col))) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = all_of(pollinator_col),
                values_from = present,
                values_fill = 0) %>%
    tibble::column_to_rownames(plant_col) %>%
    as.matrix()
  
  mat <- mat > 0
  
  # Filter to plants and pollinators of interest
  mat[rownames(mat) %in% plants,
      colnames(mat) %in% pollinators,
      drop = FALSE]
}

#### Complimentary mix wrappers ####
# Simple wrapper - find complementary plants for a single region.
# Takes pre-filtered interactions for a region, returns best plant combination.
find_complementary_plants <- function(region_data,
                                      interactions_df,
                                      totals_df,
                                      target_group = "all",
                                      n_plants = 10,
                                      n_generations = 100,
                                      population_size = 200,
                                      mutation_rate = 0.1) {
  
  region_interactions <- interactions_df %>%
    filter(JEPCODE == region_data$JEPCODE[1])
  
  if (target_group != "all") {
    region_interactions <- region_interactions %>%
      filter(sourceTaxonType == target_group)
  }
  
  all_plants <- unique(region_interactions$targetTaxonName_harm)
  
  if (length(all_plants) < n_plants) {
    n_plants <- length(all_plants)
  }
  
  interaction_matrix <- build_interaction_matrix(
    region_interactions,
    plants = all_plants,
    pollinators = unique(region_interactions$sourceTaxonName_harm)
  )
  
  ga_result <- run_genetic_algorithm(
    interaction_matrix = interaction_matrix,
    n_plants = n_plants,
    n_generations = n_generations,
    population_size = population_size,
    mutation_rate = mutation_rate
  )
  
  return(ga_result$best_solution)
}

# Full wrapper - find top mix for a region with filtered
# plant/pollinator subsets, taxon summaries, and fitness history tracking.
find_top_mix <- function(ecoregion,
                         interaction_data_subset,
                         candidate_plants,
                         focal_pollinators,
                         checklist,
                         sp_by_region,
                         n_plants,
                         n_generations = 100,
                         population_size = 200,
                         mutation_rate = 0.1) {
  
  # Filter data for the specific ecoregion
  ecoregion_data <- interaction_data_subset %>%
    filter(JEPCODE %in% c(ecoregion)) %>%
    mutate(potentialOrConfirmed = 1)
  
  # Checklist of species in this ecoregion
  ecoregion_sp <- sp_by_region %>%
    filter(JEPCODE == ecoregion) %>%
    pull(genus_species)
  
  checklist_region <- checklist %>%
    dplyr::select(genus_species, taxon) %>%
    unique() %>%
    filter(genus_species %in% ecoregion_sp)
  
  # Filter plants and pollinators to those present in this ecoregion
  plants <- candidate_plants[candidate_plants %in% checklist_region[checklist_region$taxon == "plants", ]$genus_species] %>%
    unique()
  
  pollinators <- focal_pollinators[focal_pollinators %in% checklist_region[checklist_region$taxon %in% c("bees", "moths", "butterflies", "hoverflies"), ]$genus_species] %>%
    unique()
  
  # Build interaction matrix
  interaction_matrix_filtered <- build_interaction_matrix(
    ecoregion_data,
    plants = plants,
    pollinators = pollinators
  )
  
  # Check if we have enough plants
  if (nrow(interaction_matrix_filtered) < n_plants) {
    warning(paste("Only", nrow(interaction_matrix_filtered),
                  "plants available in ecoregion", ecoregion))
    n_plants <- nrow(interaction_matrix_filtered)
  }
  
  message(paste("Using genetic algorithm with", nrow(interaction_matrix_filtered), "candidate plants..."))
  
  # Run GA
  ga_result <- run_genetic_algorithm(
    interaction_matrix = interaction_matrix_filtered,
    n_plants = n_plants,
    n_generations = n_generations,
    population_size = population_size,
    mutation_rate = mutation_rate
  )
  
  best_solution <- ga_result$best_solution
  best_fitness  <- ga_result$best_fitness
  best_solutions <- ga_result$best_solutions
  fitness_history <- ga_result$fitness_history
  
  # Calculate covered pollinators
  covered_pollinators <- apply(
    interaction_matrix_filtered[best_solution, , drop = FALSE], 2, any
  )
  
  # Taxon group summary
  all_pollinators_in_region <- checklist_region %>%
    filter(taxon %in% c("bees", "moths", "butterflies", "hoverflies"))
  
  supported_pollinator_names <- names(covered_pollinators)[covered_pollinators]
  
  taxon_totals <- all_pollinators_in_region %>%
    group_by(taxon) %>%
    summarise(total_in_region = n(), .groups = "drop")
  
  taxon_supported <- all_pollinators_in_region %>%
    filter(genus_species %in% supported_pollinator_names) %>%
    group_by(taxon) %>%
    summarise(supported_by_mix = n(), .groups = "drop")
  
  taxon_summary <- taxon_totals %>%
    left_join(taxon_supported, by = "taxon") %>%
    mutate(
      supported_by_mix = replace_na(supported_by_mix, 0),
      percent_supported = round(100 * supported_by_mix / total_in_region, 1),
      ecoregion = ecoregion
    ) %>%
    dplyr::select(ecoregion, taxon, total_in_region, supported_by_mix, percent_supported)
  
  message(paste("\n=== Taxon Group Summary for", ecoregion, "==="))
  print(taxon_summary)
  
  # Results dataframe
  top_mix <- data.frame(
    ecoregion = ecoregion,
    rank = 1:n_plants,
    scientificName = best_solution,
    cumulative_pollinators = best_fitness
  )
  
  attr(top_mix, "all_best_combinations") <- best_solutions
  attr(top_mix, "n_best_combinations")   <- length(best_solutions)
  attr(top_mix, "fitness_history")       <- fitness_history
  attr(top_mix, "final_fitness")         <- best_fitness
  
  return(list(
    top_mix      = top_mix,
    taxon_summary = taxon_summary
  ))
}