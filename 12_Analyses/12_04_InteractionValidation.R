#### Load packages ####
library(dplyr)
library(stringr)
library(tidyr)
library(sf)
library(terra)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggdist)
library(units)
library(data.table)

# Note: this script loads its own data rather than sourcing 14_00_LoadData.R

#### Read in data ####
# California ecoregions
ecoregions <- read_sf("Data_Raw/Spatial/Jepson_Ecoregions/jepson.shp") %>%
  st_make_valid()

## Checklists
checklist_cleaned <- fread("Data_Clean/Species_Checklists/checklist_cleaned.csv")

# Plant checklist
plant_checklist <- checklist_cleaned[checklist_cleaned$taxon == "plants", ] %>%
  dplyr::select(genus_species) %>%
  mutate(genus = str_extract(genus_species, "[A-Z][a-z]*")) %>%
  unique()

# Pollinator checklists
pollinator_checklists <- checklist_cleaned[checklist_cleaned$taxon %in% c("bees", "butterflies", "moths", "hoverflies"), ] %>%
  dplyr::select(genus_species) %>%
  mutate(genus = str_extract(genus_species, "[A-Z][a-z]*")) %>%
  unique()

# Occurrences
occurrences <- fread("Data_Clean/Occurrences/occurrences_clean.csv") %>%
  filter(genus_species %in% checklist_cleaned$genus_species)

# Phenology
phenology <- fread("Data_Clean/Phenology/phenology_best_available.csv") %>% distinct()

## Interactions
interactions <- fread("Data_Clean/Interactions/ints_final_clean_slim.csv")

## Spatially explicit interactions in CA
interactions$decimalLatitude  <- as.numeric(interactions$decimalLatitude)
interactions$decimalLongitude <- as.numeric(interactions$decimalLongitude)

interactions_sf <- st_as_sf(
  interactions[!is.na(interactions$decimalLatitude) &
                 !is.na(interactions$decimalLongitude), ],
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326,
  remove = FALSE
)

interactions_sf_witheco <- st_join(interactions_sf, ecoregions) %>%
  st_drop_geometry()

interactions_regional <- interactions_sf_witheco %>%
  filter(!is.na(JEPCODE))

interactions_witheco <- left_join(interactions, interactions_regional) %>%
  mutate(sourceTaxonSpeciesName = word(sourceTaxonName_harm, 2),
         targetTaxonSpeciesName = word(targetTaxonName_harm, 2))

interactions_regional_visitsFlowersOf <- filter(interactions_regional, interactionTypeName == "visitsFlowersOf")
interactions_regional_hasHost         <- filter(interactions_regional, interactionTypeName == "hasHost")

spatialInteractions_unique_visitsFlowersOf <- interactions_regional_visitsFlowersOf %>%
  filter(str_count(targetTaxonName_harm, "\\S+") >= 2,
         str_count(sourceTaxonName_harm, "\\S+") >= 2,
         !is.na(JEPCODE)) %>%
  dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE, interactionTypeName) %>%
  mutate(
    sourceTaxonName_harm    = word(sourceTaxonName_harm, start = 1, end = 2),
    targetTaxonName_harm    = word(targetTaxonName_harm, start = 1, end = 2),
    sourceTaxonGenusName    = word(sourceTaxonName_harm, 1),
    targetTaxonGenusName    = word(targetTaxonName_harm, 1),
    sourceTaxonSpeciesName  = word(sourceTaxonName_harm, 2),
    targetTaxonSpeciesName  = word(targetTaxonName_harm, 2)
  ) %>%
  filter(sourceTaxonName_harm %in% checklist_cleaned$genus_species,
         targetTaxonName_harm %in% checklist_cleaned$genus_species) %>%
  unique()

#### Nectar Interaction Predictions ####
predict_visitsFlowersOf <- read.csv("Data_Clean/Interactions/predictions_supporting_visitsFlowersOf.csv")
predict_hasHost         <- read.csv("Data_Clean/Interactions/predictions_supporting_hasHost.csv")

predict_visitsFlowersOf <- predict_visitsFlowersOf %>%
  filter(method_type == "Level1_Simple") %>%
  unique()

#### Normalize Spatial Overlap ####
sdm_example <- rast(list.files("Data_Clean/SDMs/sdm_by_taxon/plants/continuous", full.names = TRUE)[1])
res_value <- res(sdm_example)[1]
crs_value <- crs(sdm_example)

ecoregions_union_full <- ecoregions %>%
  st_union() %>%
  st_make_valid()

ecoregions_vect_full <- vect(as(ecoregions_union_full, "Spatial"))

ones_rast <- rast(ext(ecoregions_vect_full), resolution = res_value, crs = crs_value)
values(ones_rast) <- 1
ones_rast <- mask(ones_rast, ecoregions_vect_full)

ecoregions_union <- ecoregions %>%
  group_by(JEPCODE) %>%
  summarize(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

unique_jepcodes <- ecoregions_union$JEPCODE

pixel_counts <- data.frame(JEPCODE = unique_jepcodes, n_pixels = NA_integer_)

for (i in seq_along(unique_jepcodes)) {
  eco_single <- ecoregions_union %>% filter(JEPCODE == unique_jepcodes[i])
  eco_vect   <- vect(as(eco_single, "Spatial"))
  masked     <- terra::mask(terra::crop(ones_rast, eco_vect), eco_vect)
  pixel_counts$n_pixels[i] <- sum(!is.na(values(masked)))
}

predict_visitsFlowersOf <- predict_visitsFlowersOf %>%
  left_join(pixel_counts, by = "JEPCODE") %>%
  mutate(spatial_overlap_normalized = spatial_overlap / n_pixels)

predict_hasHost <- predict_hasHost %>%
  left_join(pixel_counts, by = "JEPCODE") %>%
  mutate(spatial_overlap_normalized = spatial_overlap / n_pixels)

#### Normalize temporal overlap ####
predict_visitsFlowersOf <- predict_visitsFlowersOf %>%
  mutate(temporal_overlap_normalized = temporal_overlap / plant_duration)

#### Predict flower visitation interactions ####
predictions_visitsFlowersOf <- predict_visitsFlowersOf %>%
  mutate(interactionScore = case_when(
    spatial_overlap_normalized == 0 ~ 0,
    temporal_overlap_normalized == 0 ~ 0,
    TRUE ~ spatial_overlap_normalized * temporal_overlap_normalized
  ))

predictions_visitsFlowersOf_NA <- predictions_visitsFlowersOf[is.na(predictions_visitsFlowersOf$interactionScore), ]
predictions_visitsFlowersOf    <- predictions_visitsFlowersOf[!is.na(predictions_visitsFlowersOf$interactionScore), ]

predictions_visitsFlowersOf <- predictions_visitsFlowersOf %>%
  mutate(targetTaxonGenusName = sapply(strsplit(targetTaxonName_harm, " ", fixed = TRUE), `[`, 1))

predictions_hasHost <- predict_hasHost %>%
  mutate(interactionScore = case_when(
    spatial_overlap_normalized == 0 ~ 0,
    TRUE ~ spatial_overlap_normalized
  ))

predictions_hasHost_NA <- predictions_hasHost[is.na(predictions_hasHost$interactionScore), ]
predictions_hasHost    <- predictions_hasHost[!is.na(predictions_hasHost$interactionScore), ]

predictions_hasHost <- predictions_hasHost %>%
  mutate(targetTaxonGenusName = sapply(strsplit(targetTaxonName_harm, " ", fixed = TRUE), `[`, 1))

#### Leave out validation ####
all_pollinators <- pollinator_checklists$genus_species
all_plants      <- plant_checklist$genus_species
all_JEPCODE     <- unique(ecoregions$JEPCODE)

predictions_visitsFlowersOf$targetTaxonGenusName <- sub(" .*", "", predictions_visitsFlowersOf$targetTaxonName_harm)

# Validation function with multiple null runs
run_validation_flowerVisitation <- function(interactions_witheco,
                                            predictions_visitsFlowersOf,
                                            spatialInteractions_unique_visitsFlowersOf,
                                            plant_checklist,
                                            pollinator_checklists,
                                            ecoregions,
                                            validation_prop = 0.2,
                                            n_null_runs = 100,
                                            seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Define potential interaction pairs (pol species - plant genus)
  potential_pairs <- predictions_visitsFlowersOf %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonGenusName) %>%
    distinct()
  
  # 2. Randomly select interactions to remove (validation set)
  spatialInteractions_toRemove <- spatialInteractions_unique_visitsFlowersOf %>%
    slice_sample(prop = validation_prop)
  
  # 3. Remove validation set from training data
  interactions_filteredValidation <- interactions_witheco %>%
    filter(interactionTypeName == "visitsFlowersOf") %>%
    anti_join(spatialInteractions_toRemove, 
              by = c("sourceTaxonGenusName", "targetTaxonGenusName", 
                     "sourceTaxonSpeciesName", "targetTaxonSpeciesName", 
                     "JEPCODE"))
  
  # 4. Identify which potential pairs are still allowed after removal
  ints_pol_plantgenera_filteredValidation <- interactions_filteredValidation %>%
    dplyr::select(targetTaxonGenusName, sourceTaxonName_harm) %>%
    filter(str_count(sourceTaxonName_harm, "\\S+") >= 2) %>%
    mutate(sourceTaxonName_harm = word(sourceTaxonName_harm, 1, 2))
  
  ints_pol_plantgenera_checklist_filteredValidation <- ints_pol_plantgenera_filteredValidation %>%
    filter(sourceTaxonName_harm %in% pollinator_checklists$genus_species &
             targetTaxonGenusName %in% plant_checklist$genus) %>%
    unique()
  
  potential_pairs_validation <- semi_join(potential_pairs, 
                                          ints_pol_plantgenera_checklist_filteredValidation,
                                          by = c("sourceTaxonName_harm", "targetTaxonGenusName"))
  
  # 5. Mark which predictions are still "allowed" (validationPotentialLink)
  predictions_marked <- predictions_visitsFlowersOf %>%
    mutate(validationPotentialLink = ifelse(
      paste(sourceTaxonName_harm, targetTaxonGenusName) %in% 
        paste(potential_pairs_validation$sourceTaxonName_harm, 
              potential_pairs_validation$targetTaxonGenusName),
      1, 0
    ))
  
  # 6. Add observed interactions from filtered training data
  int_observed <- interactions_filteredValidation %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
    mutate(observedInteraction = 1) %>%
    distinct()
  
  predictions_with_obs <- predictions_marked %>%
    left_join(int_observed, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    mutate(observedInteraction = replace_na(observedInteraction, 0))
  
  # 7. Calculate threshold from confirmed interactions
  confirmed_probs <- predictions_with_obs %>%
    filter(observedInteraction == 1, validationPotentialLink == 1) %>%
    pull(interactionScore)
  
  lower_threshold <- quantile(confirmed_probs, 0.25, na.rm = TRUE)
  
  # 8. Filter for predicted interactions (above threshold OR observed, AND allowed)
  filtered_interactions <- predictions_with_obs %>%
    filter((interactionScore >= lower_threshold | observedInteraction == 1),
           validationPotentialLink == 1)
  
  # 9. Calculate actual recovery rate
  matches_fromRemoved <- spatialInteractions_toRemove %>%
    semi_join(filtered_interactions, 
              by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE"))
  
  actual_recovery_prop <- nrow(matches_fromRemoved) / nrow(spatialInteractions_toRemove)
  
  # 10. Count reasons for non-recovery
  n_removed_forbidden <- spatialInteractions_toRemove %>%
    anti_join(potential_pairs_validation,
              by = c("sourceTaxonName_harm")) %>%
    nrow()
  
  n_removed_below_threshold <- spatialInteractions_toRemove %>%
    semi_join(potential_pairs_validation, 
              by = c("sourceTaxonName_harm")) %>%
    anti_join(filtered_interactions,
              by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    nrow()
  
  # 11. Null models - run multiple times
  n_draw <- filtered_interactions %>%
    anti_join(int_observed, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE")) %>%
    nrow()
  
  all_pollinators <- pollinator_checklists$genus_species
  all_plants <- plant_checklist$genus_species
  all_JEPCODE <- unique(ecoregions$JEPCODE)
  
  # Initialize vectors
  null_1_props <- numeric(n_null_runs)
  null_2_props <- numeric(n_null_runs)
  
  # For constrained validation (Null 2), what possible pairs?
  potential_for_null2 <- predictions_with_obs %>%
    filter(validationPotentialLink == 1,
           observedInteraction == 0) %>%
    dplyr::select(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %>%
    distinct()
  
  for (null_run in 1:n_null_runs) {
    ## Null model 1 - completely random interactions
    null_1 <- expand.grid(
      sourceTaxonName_harm = sample(all_pollinators, min(length(all_pollinators), 500)),
      targetTaxonName_harm = sample(all_plants, min(length(all_plants), 500)),
      JEPCODE = all_JEPCODE,
      stringsAsFactors = FALSE
    ) %>%
      slice_sample(n = n_draw, replace = FALSE)
    
    null_1_overlap <- spatialInteractions_toRemove %>%
      semi_join(null_1, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE"))
    
    null_1_props[null_run] <- nrow(null_1_overlap) / nrow(spatialInteractions_toRemove)
    
    ## Null model 2 - random from allowed potential interactions
    if (nrow(potential_for_null2) >= n_draw) {
      null_2 <- potential_for_null2 %>%
        slice_sample(n = n_draw, replace = FALSE)
      
      null_2_overlap <- spatialInteractions_toRemove %>%
        semi_join(null_2, by = c("sourceTaxonName_harm", "targetTaxonName_harm", "JEPCODE"))
      
      null_2_props[null_run] <- nrow(null_2_overlap) / nrow(spatialInteractions_toRemove)
    } else {
      null_2_props[null_run] <- NA
    }
  }
  
  # 12. Calculate p-values
  null_1_beats_actual <- mean(null_1_props >= actual_recovery_prop)
  null_2_beats_actual <- mean(null_2_props >= actual_recovery_prop, na.rm = TRUE)
  
  # 13. Return results
  list(
    summary = data.frame(
      validation_prop = validation_prop,
      n_removed = nrow(spatialInteractions_toRemove),
      n_predicted = nrow(filtered_interactions),
      n_observed_in_training = nrow(int_observed),
      threshold = lower_threshold,
      actual_overlap_n = nrow(matches_fromRemoved),
      actual_recovery_prop = actual_recovery_prop,
      n_removed_forbidden = n_removed_forbidden,
      n_removed_below_threshold = n_removed_below_threshold,
      null1_mean_prop = mean(null_1_props),
      null1_sd_prop = sd(null_1_props),
      null1_beats_actual_prop = null_1_beats_actual,
      null2_mean_prop = mean(null_2_props, na.rm = TRUE),
      null2_sd_prop = sd(null_2_props, na.rm = TRUE),
      null2_beats_actual_prop = null_2_beats_actual
    ),
    null_1_props = null_1_props,
    null_2_props = null_2_props,
    removed_interactions = spatialInteractions_toRemove,
    recovered_interactions = matches_fromRemoved
  )
}

# wrapper function
run_multiple_validation_flowerVisitation <- function(interactions_witheco,
                                                     predictions_visitsFlowersOf,
                                                     spatialInteractions_unique_visitsFlowersOf,
                                                     plant_checklist,
                                                     pollinator_checklists,
                                                     ecoregions,
                                                     n_runs = 100,
                                                     n_null_runs = 100,
                                                     validation_prop = 0.2) {
  
  # Initialize results list
  results_list <- list()
  all_null_1_props <- list()
  all_null_2_props <- list()
  
  # Run validation in loop with progress
  for (i in 1:n_runs) {
    cat("Completed", i, "of", n_runs, "runs\n")
    
    result <- run_validation_flowerVisitation(
      interactions_witheco = interactions_witheco,
      predictions_visitsFlowersOf = predictions_visitsFlowersOf,
      spatialInteractions_unique_visitsFlowersOf = spatialInteractions_unique_visitsFlowersOf,
      plant_checklist = plant_checklist,
      pollinator_checklists = pollinator_checklists,
      ecoregions = ecoregions,
      validation_prop = validation_prop,
      n_null_runs = n_null_runs,
      seed = i
    )
    
    results_list[[i]] <- result$summary
    all_null_1_props[[i]] <- result$null_1_props
    all_null_2_props[[i]] <- result$null_2_props
  }
  
  # Combine results
  list(
    summary_df = do.call(rbind, results_list),
    null_1_distributions = all_null_1_props,
    null_2_distributions = all_null_2_props
  )
}

# Run the validation
validation_results_flowerVisitation <- run_multiple_validation_flowerVisitation(
  interactions_witheco = interactions_witheco,
  predictions_visitsFlowersOf = predictions_visitsFlowersOf,
  spatialInteractions_unique_visitsFlowersOf = spatialInteractions_unique_visitsFlowersOf,
  plant_checklist = plant_checklist,
  pollinator_checklists = pollinator_checklists,
  ecoregions = ecoregions,
  n_runs = 25,
  n_null_runs = 100,
  validation_prop = 0.2
)

# View summary statistics across all runs
# summary(validation_results_flowerVisitation$summary_df)

saveRDS(validation_results_flowerVisitation, "Data_Clean/Validation/interaction_validation_results.rds")

validation_results_flowerVisitation <- readRDS("Data_Clean/Validation/interaction_validation_results.rds")

### Plot: Interaction plausibility score histograms ###
predictions_visitsFlowersOf <- predictions_visitsFlowersOf %>%
  mutate(observed = paste(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %in%
           paste(spatialInteractions_unique_visitsFlowersOf$sourceTaxonName_harm,
                 spatialInteractions_unique_visitsFlowersOf$targetTaxonName_harm,
                 spatialInteractions_unique_visitsFlowersOf$JEPCODE))



# Split data
obs_true  <- predictions_visitsFlowersOf %>% filter(observed == TRUE)
obs_false <- predictions_visitsFlowersOf %>% filter(observed == FALSE)

scale_factor <- max(table(cut(obs_false$interactionScore, breaks = 30))) /
  max(table(cut(obs_true$interactionScore,  breaks = 30)))

hist_scores <- ggplot() +
  geom_histogram(data = obs_false,
                 aes(x = interactionScore),
                 bins = 30,
                 fill = "#FFCC00", color = "#331E0F") +
  geom_histogram(data = obs_true,
                 aes(x = interactionScore, y = after_stat(count) * scale_factor),
                 bins = 30,
                 fill = "#BAFF76", color = "#3B7004", alpha = 0.5) +
  geom_vline(xintercept = quantile(obs_true$interactionScore, 0.25),
              linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_y_continuous(
    name = "Density of unobserved\ninteractions (yellow)",
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Density of observed spatially-\nexplicit interactions (green)")
  ) +
  scale_x_continuous(name = "Interaction Score") +
  theme_cowplot()

## Has Host Interaction Scores Plot
### Plot: Interaction plausibility score histograms ###
predictions_hasHost <- predictions_hasHost %>%
  mutate(observed = paste(sourceTaxonName_harm, targetTaxonName_harm, JEPCODE) %in%
           paste(spatialInteractions_unique_visitsFlowersOf$sourceTaxonName_harm,
                 spatialInteractions_unique_visitsFlowersOf$targetTaxonName_harm,
                 spatialInteractions_unique_visitsFlowersOf$JEPCODE))

# Split data
obs_true  <- predictions_hasHost %>% filter(observed == TRUE)
obs_false <- predictions_hasHost %>% filter(observed == FALSE)

scale_factor <- max(table(cut(obs_false$interactionScore, breaks = 30))) /
  max(table(cut(obs_true$interactionScore,  breaks = 30)))

hist_scores_hasHost <- ggplot() +
  geom_histogram(data = obs_false,
                 aes(x = interactionScore),
                 bins = 30,
                 fill = "#FFCC00", color = "#331E0F") +
  geom_histogram(data = obs_true,
                 aes(x = interactionScore, y = after_stat(count) * scale_factor),
                 bins = 30,
                 fill = "#BAFF76", color = "#3B7004", alpha = 0.5) +
  geom_vline(xintercept = quantile(obs_true$interactionScore, 0.25),
             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_y_continuous(
    name = "Density of unobserved\ninteractions (yellow)",
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Density of observed spatially-\nexplicit interactions (green)")
  ) +
  scale_x_continuous(name = "Interaction Score") +
  theme_cowplot()

combined_hist <- (hist_scores | hist_scores_hasHost) +
  plot_annotation(tag_levels = 'A') &
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))


ggsave("Figures/Supplementary/prediction_interactionScore_histogram.png", 
       combined_hist, 
       width = 200, 
       height = 100, 
       units = "mm", 
       dpi = 2000)

# corresponding percentile for unobserved
mean(obs_false$interactionScore <= quantile(obs_true$interactionScore, 0.25))


### Plot: Recall rate comparison null models vs NECTAR
# Prepare data - expand null distributions
plot_data <- bind_rows(
  # Null 1 distributions
  tibble(
    fold = rep(1:25, each = 100),
    model = "Null 1: Completely Random",
    recovery_prop = unlist(validation_results_flowerVisitation$null_1_distributions)
  ),
  # Null 2 distributions
  tibble(
    fold = rep(1:25, each = 100),
    model = "Null 2: Phylogenetically Constrained",
    recovery_prop = unlist(validation_results_flowerVisitation$null_2_distributions)
  )
) %>%
  mutate(fold = factor(fold, levels = 25:1))  # Reverse for top-to-bottom

# Actual recovery as points
actual_data <- validation_results_flowerVisitation$summary_df %>%
  mutate(
    fold = factor(row_number(), levels = 25:1),
    model = "Actual"
  ) %>%
  dplyr::select(fold, model, recovery_prop = actual_recovery_prop)

validation_recall <- ggplot() +
  # show intervals of null distributions
  stat_interval(data = plot_data,
                aes(y = fold, x = recovery_prop, color = model),
                .width = c(.50, .95),
                position = position_dodge(width = 0.4),
                linewidth = 2) +
  # add points for actual recovery (now with shape aesthetic mapped)
  geom_point(data = actual_data,
             aes(y = fold, x = recovery_prop, shape = model),
             size = 3, color = "#FFCC00") +
  geom_vline(xintercept = 0, linetype = "solid", alpha = 0.3) +
  scale_color_manual(
    values = c("Null 1: Completely Random" = "gray70", 
               "Null 2: Phylogenetically Constrained" = "gray40"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Actual" = 18),
    name = NULL,
    labels = c("NECTAR")
  ) +
  scale_y_discrete(labels = 1:25) +
  labs(
    x = "Recall Rate (Proportion of withheld interactions recovered)",
    y = "Validation Fold"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90")
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

ggsave("Figures/Supplementary/validation_recall.png", 
       validation_recall, 
       width = 200, 
       height = 200, 
       units = "mm", 
       dpi = 2000)

## Plot: Summary of why non-recovery of interaction
reason_summary <- validation_results_flowerVisitation$summary_df %>%
  summarise(
    `Recovered` = sum(actual_overlap_n),
    `Forbidden interaction` = sum(n_removed_forbidden),
    `Below threshold` = sum(n_removed_below_threshold)
  ) %>%
  pivot_longer(everything(), names_to = "reason", values_to = "n_interactions") %>%
  mutate(
    pct = round(n_interactions / sum(n_interactions) * 100, 1),
    reason = factor(reason, levels = c("Recovered", "Below threshold", "Forbidden interaction"))
  )

# Create pie chart
validation_outcome <- ggplot(reason_summary, aes(x = "", y = n_interactions, fill = reason)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
  coord_polar("y") +
  scale_fill_manual(
    values = c("Recovered" = "#FFCC00", 
               "Below threshold" = "gray60", 
               "Forbidden interaction" = "gray80")
  ) +
  theme_void() +
  theme(legend.position = "right") +
  labs(
    fill = "Outcome of\nWithheld Interaction",
  ) +
  # add percentage labels
  geom_text(aes(label = paste0(pct, "%")),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold")

ggsave("Figures/Supplementary/validation_outcomeBreakdown.png", 
       validation_outcome, 
       width = 200, 
       height = 200, 
       units = "mm", 
       dpi = 2000)
