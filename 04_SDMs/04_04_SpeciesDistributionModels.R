#### Array job information ####

#Right now the data is split into 10 roughly equal sections in "/input/occ_split" and in this script, each portion is called by setting the integer f (which is used to define the occ_path below). So if f = 1, it will use the first portion of data, etc. Can run the whole thing by running: 

# for (f in 1:11) {
#   source("path_to_this_script.R")
# }

#### SDM Script ####

# SDMs for native plants and insects in California
# For: Morpho/CNPS/Calscape
# By: Chris Cosma, Conservation Biology Institute
# Revised 2025-11-07 

#### Setup ####
set.seed(123)

# Command line args for array job
args <- commandArgs(trailingOnly = TRUE)
f <- as.numeric(args[1])

#Using covsel package for automated predictor selection: https://www.sciencedirect.com/science/article/pii/S1574954123001097
#needs to be installed via GitHub
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

if (!requireNamespace("covsel", quietly = TRUE)) {
  remotes::install_github("antadde/covsel", ref = "v1.0.0", upgrade = FALSE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(sf)
  library(raster) 
  library(flexsdm)
  library(biomod2)
  library(ecospat)
  library(blockCV)
  library(covsel) 
  library(dismo)  
  library(adehabitatHR)  
  library(mda)    
  library(earth)
  library(progress)
  library(rnaturalearth)
})

#disable all plotting during this session
options(device = function(...) grDevices::pdf(NULL))

# Helpers
`%||%` <- function(a, b) if (is.null(a)) b else a

# CRS
ca_albers <- "EPSG:3310"

# Paths

#set working directory
# setwd("~/Desktop")

# Define root working directory 
root_dir <- getwd() 

# Define and create sub directories

#this should already be created and contain all the data
input_dir = file.path(root_dir, "Data_Clean/SDMs/SDMs_Inputs")

#these can be made here
intermediate_dir = file.path(root_dir, "Temp")
dir.create(intermediate_dir, showWarnings = FALSE, recursive = TRUE)

output_dir = file.path(root_dir, "Data_Clean/SDMs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

biomod_output_dir <- file.path(intermediate_dir, "biomod2_models")
dir.create(biomod_output_dir, showWarnings = FALSE, recursive = TRUE)

# Create taxon-specific output directories
taxon_output_base <- file.path(output_dir, "sdm_by_taxon")
dir.create(taxon_output_base, showWarnings = FALSE, recursive = TRUE)

#define file paths
tif_path = file.path(input_dir, "env_stack_final.tif")
eco_path = file.path(input_dir, "jepson_epa_ecoregions_combined.gpkg")

#Define occurrence data portion path (for array job)
occ_path = file.path(input_dir, paste0("occ_split/sdm_data_", f, ".csv"))

#### Load data ####
occ_data <- read.csv(occ_path)

# Get unique taxa for output organization
unique_taxa <- unique(occ_data$taxon)
message(sprintf("Found %d unique taxa: %s", length(unique_taxa), paste(unique_taxa, collapse = ", ")))

# Create output directories for each taxon
for (taxon in unique_taxa) {
  taxon_dir <- file.path(taxon_output_base, gsub("[^[:alnum:]_-]", "_", taxon))
  dir.create(file.path(taxon_dir, "continuous"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(taxon_dir, "binary"), showWarnings = FALSE, recursive = TRUE)
}

#load combined ecoregion file (from script "joined_exoregions.R")
jepson_proj = st_read(eco_path, quiet = TRUE) |> 
  st_make_valid() |>
  st_transform(ca_albers)

ca_boundary <- readRDS(file.path(input_dir,"/gadm41_USA_1_pk.rds"))
ca_boundary <- ca_boundary[ca_boundary$NAME_1 == "California", ]
ca_boundary_proj <- st_transform(st_as_sf(ca_boundary), ca_albers)

# Environmental stack
env_stack <- terra::rast(tif_path)

# Keep only valid layers
valid_layers <- which(terra::global(env_stack, "max", na.rm = TRUE)[, 1] != -Inf)
env_stack <- env_stack[[valid_layers]]

# CA mask/template
ca_bnd_wgs <- st_transform(ca_boundary_proj, 4326) |> terra::vect()
template   <- terra::mask(terra::crop(env_stack[[1]], ca_bnd_wgs), ca_bnd_wgs)

#### Categorize Predictors ####

map_layer_category <- function(nm) {
  nm_l <- tolower(nm)
  dplyr::case_when(
    # All climate variables (temp, precip, water balance, radiation)
    grepl("bio|cwd|aet|pet|vpd|srad|solar|gdd|frost|terraclim", nm_l) ~ "climate",
    
    # All topographic variables (elevation, slope, aspect, terrain complexity)
    grepl("elev|slope|aspect|tri|roughness|rugged|terrain", nm_l) ~ "topography",
    
    # Vegetation state & productivity (NDVI, LAI, land cover)
    grepl("ndvi|evi|lai|leaf|^lc_|land.*cover|nlcd.*mode", nm_l) ~ "vegetation",
    
    # All soil properties (texture, chemistry, organic matter)
    grepl("^soil", nm_l) ~ "soil",
    
    # Hydrologic features (distance to water, streams, coast)
    grepl("dist.*coast|dist.*stream|dist.*water", nm_l) ~ "hydrology",
    
    # Human impact & disturbance (fragmentation, fire, development)
    grepl("human.*footprint|landscape.*intactness|impervious|pop2020|fire.*freq|time.*since.*fire|burn", nm_l) ~ "disturbance",
    
    # Catch-all
    TRUE ~ "other"
  )
}

# Apply categorization
layer_names <- names(env_stack)
layer_categories <- purrr::map_chr(layer_names, map_layer_category)

# Summary
cat("\n=== Predictor Category Summary ===\n")
print(table(layer_categories))

# Create reference table
category_reference <- data.frame(
  layer = layer_names,
  category = layer_categories,
  stringsAsFactors = FALSE
) %>%
  arrange(category, layer)


cat("\nCategory distribution:\n")
print(
  category_reference %>%
    count(category, name = "n_predictors") %>%
    arrange(desc(n_predictors))
)

# key categories by layer name so we can align reliably later
names(layer_categories) <- layer_names

message(sprintf("Loaded %d environmental predictors across %d categories", 
                length(layer_names), length(unique(layer_categories))))

#### One-per-cell thinning ####
one_per_cell <- function(df, template) {
  if (!nrow(df)) return(df[0, ])
  cells <- terra::cellFromXY(template, as.matrix(df[, c("longitude", "latitude")]))
  keep  <- !duplicated(cells)
  df[keep, , drop = FALSE]
}

#### Target Group Background Sampling ####
generate_target_group_background <- function(env_masked, accessible_vect, target_species, occ_data, n_bg_target = 5000L) {
  # 1) Same taxon, other species
  target_taxon <- unique(occ_data$taxon[
    occ_data$species == target_species | occ_data$species_full == target_species
  ])
  tg <- occ_data %>%
    filter(taxon == target_taxon, 
           !(species == target_species | species_full == target_species)) %>%
    dplyr::select(longitude, latitude) %>%
    distinct()
  
  # 2) Fallback to all others if too few
  if (nrow(tg) < n_bg_target) {
    tg <- occ_data %>%
      filter(species != target_species) %>%
      dplyr::select(longitude, latitude) %>%
      distinct()
  }
  
  # 3) Restrict to accessible area
  if (!inherits(accessible_vect, "SpatVector")) accessible_vect <- terra::vect(accessible_vect)
  tg_sf <- sf::st_as_sf(tg, coords = c("longitude", "latitude"), crs = 4326)
  inside <- sf::st_within(tg_sf, sf::st_as_sf(accessible_vect), sparse = FALSE)[, 1]
  tg <- tg[inside, , drop = FALSE]
  
  # 4) Keep only locations with valid env data
  tg_env <- tryCatch({
    terra::extract(env_masked, as.matrix(tg[, c("longitude", "latitude")]))
  }, error = function(e) NULL)
  
  if (is.null(tg_env) || nrow(tg_env) == 0) {
    tg <- tg[0, , drop = FALSE]
  } else {
    keep <- stats::complete.cases(tg_env)
    tg <- tg[keep, , drop = FALSE]
  }
  
  # 5) One-per-cell thinning
  if (nrow(tg) > 0) {
    cells <- terra::cellFromXY(env_masked[[1]], as.matrix(tg))
    tg <- tg[!duplicated(cells), , drop = FALSE]
  }
  
  # 6) Random pseudo-absence fallback if empty
  n_av <- nrow(tg)
  if (n_av == 0) {
    warning("No valid target-group background; using random pseudo-absences inside accessible area.")
    rand_bg <- terra::spatSample(env_masked[[1]], n = n_bg_target, method = "random", na.rm = TRUE, xy = TRUE)
    rand_bg <- as.data.frame(rand_bg)[, c("x", "y")]
    colnames(rand_bg) <- c("longitude", "latitude")
    return(rand_bg)
  }
  
  # 7) Sample to target size
  tg[sample(n_av, min(n_bg_target, n_av), replace = n_av < n_bg_target), , drop = FALSE]
}

#### Threshold helpers ####

#For now just using OR10
THRESHOLD_METHOD <- "OR10"

# calculate_f_threshold <- function(obs, pred) {
#   keep <- !(is.na(obs) | is.na(pred)); if (sum(keep) < 2) return(0.5)
#   obs <- obs[keep]; pred <- pred[keep]; if (length(unique(obs)) < 2) return(0.5)
#   thr <- seq(0, 1, by = 0.01); fval <- numeric(length(thr))
#   for (i in seq_along(thr)) {
#     pb <- as.integer(pred >= thr[i])
#     tp <- sum(obs == 1 & pb == 1); fp <- sum(obs == 0 & pb == 1); fn <- sum(obs == 1 & pb == 0)
#     precision <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
#     recall    <- ifelse((tp + fn) > 0, tp / (tp + fn), 0)
#     fval[i] <- ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), 0)
#   }
#   out <- thr[which.max(fval)]; if (!is.finite(out)) 0.5 else out
# }
# 
# calculate_mcc_threshold <- function(obs, pred) {
#   keep <- !(is.na(obs) | is.na(pred)); if (sum(keep) < 2) return(0.5)
#   obs <- obs[keep]; pred <- pred[keep]; if (length(unique(obs)) < 2) return(0.5)
#   thr <- seq(0, 1, by = 0.01); mcc_val <- numeric(length(thr))
#   for (i in seq_along(thr)) {
#     pb <- as.integer(pred >= thr[i])
#     tp <- sum(obs == 1 & pb == 1); fp <- sum(obs == 0 & pb == 1); tn <- sum(obs == 0 & pb == 0); fn <- sum(obs == 1 & pb == 0)
#     denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
#     mcc_val[i] <- ifelse(denom > 0, (tp * tn - fp * fn) / denom, 0)
#   }
#   out <- thr[which.max(mcc_val)]; if (!is.finite(out)) 0.5 else out
# }

calculate_OR10_threshold <- function(obs, pred, min_thr = 0.01) {
  keep <- !(is.na(obs) | is.na(pred))
  if (sum(keep) < 2) return(0.5)
  obs <- obs[keep]; pred <- pred[keep]
  p1 <- obs == 1
  if (!any(p1)) return(0.5)
  pred_pres <- pred[p1]
  thr <- quantile(pred_pres, probs = 0.10, na.rm = TRUE)
  max(thr, min_thr)
}


#### Core function for building SDMs ####
sdm_for_species <- function(species_name, occ_data, env_stack, layer_categories, 
                            category_reference, jepson_proj, ca_boundary_proj,
                            biomod_output_dir) {
  
  require(flexsdm)
  require(biomod2)
  require(ecospat)
  require(blockCV)
  require(covsel)
  require(raster)  
  require(mda)    
  require(earth)
  
  # SET WORKING DIRECTORY TO BIOMOD OUTPUT DIR
  setwd(biomod_output_dir)
  
  # Retrieve taxon for this species
  species_taxon <- unique(occ_data$taxon[
    occ_data$species == species_name | occ_data$species_full == species_name
  ])
  
  # Step 1: Filter species occurrences
  sp_data <- occ_data %>%
    filter(species == species_name | species_full == species_name) %>%
    mutate(species = species_name) %>%        # ensure the modeled name is consistent if the subsepcies in species_full was selected
    dplyr::select(species, longitude, latitude)
  
  #define # of original occurrences
  n_occ_orig <- nrow(sp_data)
  
  # Step 2: One-per-cell thinning (ORIGINAL APPROACH RESTORED)
  sp_thin <- one_per_cell(sp_data, env_stack[[1]])
  n_occ_thin <- nrow(sp_thin)
  
  message(sprintf("  After one-per-cell thinning: %d occurrences", n_occ_thin))
  
  #define the raw occurrence number, which defines whether the species is ultra-rare, rare, or common. Criteria: If sp_thin >= 6, use n_occ_orig, else use n_occ_thin. Because it there are less than 6 unique environmental conditions for the model to run, its not going to run well.
  n_occ_raw <- ifelse(n_occ_thin >= 6, n_occ_orig, n_occ_thin)
  message(sprintf("\n=== Processing %s (n=%d occurrences) ===", species_name, n_occ_raw))
  
  # Step 3: Environmental outlier filtering 
  # Determine if species is rare or ultra-rare based on count above
  is_rare_or_ultra <- n_occ_raw <= 24  # Rare: 7-24, Ultra-rare: <7
  
  if (!is_rare_or_ultra) {
    # Only do environmental filtering for regular species
    sp_env <- terra::extract(env_stack, sp_thin[, c("longitude", "latitude")])
    if (sum(stats::complete.cases(sp_env[, -1])) < 3) {
      message("  → Too few occurrences with complete environmental data, skipping environmental outlier filtering")
      sp_clean <- sp_thin
    } else {
      # NOW it only runs if there ARE enough complete cases
      sp_clean <- tryCatch({
        flexsdm::occfilt_env(
          data = sp_thin, x = "longitude", y = "latitude", id = "species",
          env_layer = env_stack, nbins = 5
        )
      }, error = function(e) {
        message("  Environmental filtering failed, using thinned occurrences")
        sp_thin
      })
    }
    n_occ <- nrow(sp_clean)
  } else {
    # Skip environmental filtering for rare/ultra-rare species
    message(sprintf("  Species is %s (n=%d), skipping environmental outlier filtering", 
                    ifelse(n_occ_raw < 7, "ultra-rare", "rare"), n_occ_raw))
    sp_clean <- sp_thin
    n_occ <- nrow(sp_clean)
  }
  
  message(sprintf("  Final occurrence count: %d", n_occ))
  
  # Step 4: Create accessible area based on combined ecoregions (CA + buffer)
  sp_points <- st_as_sf(sp_clean, coords = c("longitude", "latitude"), crs = 4326)
  sp_points_proj <- st_transform(sp_points, ca_albers)
  
  # Identify unique ecoregions containing occurrences (includes both Jepson and EPA)
  eco_intersect <- st_intersection(sp_points_proj, jepson_proj)
  unique_ecoregions <- unique(eco_intersect$ecoregion_code)
  
  # Create FULL accessible area for model training (includes CA + buffer ecoregions)
  ecoregion_mask <- jepson_proj %>% 
    filter(ecoregion_code %in% unique_ecoregions) %>%
    st_union()
  
  # IMPORTANT: Create strict mask = ONLY Jepson (CA) portions of selected ecoregions
  # This ensures models train on full extent but outputs are CA-only
  strict_ecoregion_mask <- jepson_proj %>% 
    filter(ecoregion_code %in% unique_ecoregions,
           source == "Jepson") %>%  # Only CA portions
    st_union() %>%
    st_make_valid() %>%  # VALIDATE BEFORE transforming
    st_buffer(dist = 0) %>%  # CLEAN BEFORE transforming
    st_transform(4326) %>%
    st_make_valid() %>%  # VALIDATE AFTER transforming
    st_buffer(dist = 0)  # CLEAN AFTER transforming
  
  strict_ecoregion_vect <- terra::vect(strict_ecoregion_mask)
  
  # Buffer for accessible area (applied to combined CA + buffer ecoregions)
  buffer_dist <- case_when(
    n_occ < 7 ~ 25000,   # 25 km for ultra-rare
    n_occ < 30 ~ 50000,  # 50 km for rare
    TRUE ~ 100000        # 100 km for common
  )
  
  # Create buffered accessible area - includes full extent (CA + buffer)
  accessible_area <- st_buffer(ecoregion_mask, dist = buffer_dist) %>%
    st_union()
  
  accessible_area <- st_transform(accessible_area, 4326)
  accessible_vect <- terra::vect(accessible_area)
  
  # Mask environmental data to accessible area (now includes buffer zone)
  env_masked <- terra::mask(terra::crop(env_stack, accessible_vect), accessible_vect)
  
  # Step 5: Model based on sample size 
  model_type <- NA
  final_prediction <- NULL
  performance_metrics <- list()
  n_predictors_selected <- NA
  selected_method <- NA
  selected_predictors <- NULL
  
  if (n_occ < 7) {
    # Ultra-rare: NO PREDICTOR SELECTION - use all available in accessible area
    message("  Building alpha-hull model for ultra-rare species...")
    model_type <- "ultraRare_alphahull"
    
    # Use all available environmental layers (no selection)
    env_selected <- env_masked
    n_predictors_selected <- terra::nlyr(env_selected)
    selected_method <- "none"
    selected_predictors <- names(env_selected)
    
    # Check if we have enough points for hull creation
    if (n_occ < 3) {
      # For 1–2 points, use simple buffer approach
      message("    Using buffer approach for n < 3 occurrences")
      
      # Compute interpoint distance if there are two points
      buffer_dist_ultra <- ifelse(n_occ == 1, 15000, 10000)  # 15km for 1 point, 10km for 2
      
      if (n_occ == 2) {
        coords_matrix <- coordinates(as(sp_points_proj, "Spatial"))
        dist_m <- as.numeric(sp::spDists(coords_matrix)[1, 2])
        if (dist_m > 50000) {
          message(sprintf("    Two occurrences are %.1f km apart → creating separate buffers", dist_m / 1000))
          # Keep buffers separate (do NOT union)
          buf_list <- lapply(1:2, function(i) st_buffer(sp_points_proj[i, ], dist = buffer_dist_ultra))
          pts_buffer <- do.call(sf::st_union, buf_list) |> st_transform(4326)
        } else {
          message(sprintf("    Two occurrences are only %.1f km apart → merging into single buffer", dist_m / 1000))
          pts_buffer <- st_buffer(sp_points_proj, dist = buffer_dist_ultra) |> st_union() |> st_transform(4326)
        }
      } else {
        # Single occurrence → simple 15 km buffer
        pts_buffer <- st_buffer(sp_points_proj, dist = buffer_dist_ultra) |> st_transform(4326)
      }
      
      # Clip to strict ecoregions
      range_mask_strict <- st_intersection(pts_buffer, strict_ecoregion_mask)
      buf_vect <- terra::vect(range_mask_strict)
    } else {
      # For 3–6 points: decide whether to build a hull or just buffers
      coords_matrix <- coordinates(as(sp_points_proj, "Spatial"))
      pts_spatial <- sp::SpatialPoints(
        coords_matrix,
        proj4string = sp::CRS(as.character(st_crs(sp_points_proj)$proj4string))
      )
      
      # Compute pairwise distances (meters, in Albers)
      dist_mat <- as.matrix(sp::spDists(pts_spatial))
      max_dist <- max(dist_mat[upper.tri(dist_mat)], na.rm = TRUE)
      
      # Threshold: if points are too far apart, use buffers instead of hull
      # e.g., > 50 km apart ⇒ treat as separate subpopulations
      if (max_dist > 50000) {
        message(sprintf("    Points are widely separated (max %.1f km); using simple per-point buffers", max_dist / 1000))
        pts_buffer <- st_buffer(sp_points_proj, dist = 10000) |>
          st_union() |>
          st_transform(4326)
        range_mask_strict <- st_intersection(pts_buffer, strict_ecoregion_mask)
        buf_vect <- terra::vect(range_mask_strict)
      } else {
        # Points are close enough ⇒ build a LoCoH/MCP hull
        message(sprintf("    Points are spatially clustered (max %.1f km); building LoCoH hull", max_dist / 1000))
        hull_result <- tryCatch({
          k_value <- max(2, min(n_occ - 1, 5))
          adehabitatHR::LoCoH.k(pts_spatial, k = k_value)
        }, error = function(e) {
          message("    LoCoH failed, using MCP fallback")
          tryCatch({
            adehabitatHR::mcp(pts_spatial, percent = 100)
          }, error = function(e2) {
            message("    MCP also failed, using buffer fallback")
            NULL
          })
        })
        
        if (!is.null(hull_result)) {
          hull_sf <- st_as_sf(hull_result) |> st_transform(4326)
          range_mask_strict <- st_intersection(hull_sf, strict_ecoregion_mask)
          buf_vect <- terra::vect(range_mask_strict)
        } else {
          # Final fallback if hull creation failed
          pts_buffer <- st_buffer(sp_points_proj, dist = 10000) |>
            st_union() |>
            st_transform(4326)
          range_mask_strict <- st_intersection(pts_buffer, strict_ecoregion_mask)
          buf_vect <- terra::vect(range_mask_strict)
        }
      }
    }
    
    # Create prediction with distance decay
    # Extract coordinates properly from sp_points
    coords_df <- as.data.frame(st_coordinates(sp_points))
    names(coords_df) <- c("longitude", "latitude")
    occ_vect <- terra::vect(coords_df, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    
    dist_raster <- terra::distance(env_selected[[1]], occ_vect)
    
    # Normalize distance and invert
    dist_norm <- dist_raster / terra::global(dist_raster, "max", na.rm = TRUE)[1, 1]
    final_prediction <- 1 - dist_norm
    
    # Mask to alpha-hull/buffer area
    final_prediction <- terra::mask(final_prediction, buf_vect)
    
    performance_metrics <- list(
      model_type = "ultraRare_alphahull",
      n_occurrences = n_occ,
      hull_type = ifelse(n_occ < 3, paste0("buffer_", ifelse(n_occ == 1, "15km", "10km")), 
                         ifelse(exists("k_value"), paste0("LoCoH_k", k_value), "MCP")),
      n_predictors = NA
    )
    
  } else if (n_occ <= 24) {
    # Rare: ESM with ONLY COLLINEARITY FILTERING (skip embedding)
    message("  Building ESM for rare species...")
    model_type <- "ESM"
    
    # Predictor selection - COLLINEARITY FILTERING ONLY for ESM
    message("  Running collinearity filtering only (ESM optimization)...")
    
    sp_env_df <- terra::extract(env_masked, sp_clean[, c("longitude", "latitude")])
    sp_env_df <- sp_env_df[, -1]  # Remove ID column
    
    # Remove predictors with zero variance
    var_check <- apply(sp_env_df, 2, function(x) var(x, na.rm = TRUE))
    keep_vars <- names(var_check)[var_check > 0 & !is.na(var_check)]
    sp_env_df <- sp_env_df[, keep_vars, drop = FALSE]
    env_masked_subset <- env_masked[[keep_vars]]
    
    # Store original names before cleaning
    original_names <- names(env_masked_subset)
    
    # Clean column names for covsel
    sp_env_df_clean <- sp_env_df
    colnames(sp_env_df_clean) <- gsub("[^A-Za-z0-9_]", "_", colnames(sp_env_df_clean))
    colnames(sp_env_df_clean) <- make.unique(colnames(sp_env_df_clean), sep = "_")
    
    # Get categories for the layers
    layer_categories_subset <- layer_categories[original_names]
    layer_categories_clean <- layer_categories_subset
    names(layer_categories_clean) <- colnames(sp_env_df_clean)
    
    # ONLY do collinearity filtering with covsel.filter (NO embedding for rare species)
    selected_predictors <- tryCatch({
      filtered <- covsel::covsel.filter(
        covdata = sp_env_df_clean,
        pa = rep(1, nrow(sp_env_df_clean)),
        corcut = 0.7,
        categories = layer_categories_clean[colnames(sp_env_df_clean)]
      )
      
      # Get filtered variable names
      if (is.character(filtered)) {
        # Map back to original names
        idx <- which(colnames(sp_env_df_clean) %in% filtered)
        original_names[idx]
      } else if (is.data.frame(filtered)) {
        idx <- which(colnames(sp_env_df_clean) %in% colnames(filtered))
        original_names[idx]
      } else {
        NULL
      }
    }, error = function(e) {
      message(sprintf("    Collinearity filtering failed: %s, using top 6 by variance", e$message))
      NULL
    })
    
    # Simple fallback if covsel.filter fails or returns too few
    if (is.null(selected_predictors) || length(selected_predictors) < 3) {
      message("    Using top 6 predictors by variance (fallback)")
      var_scores <- apply(sp_env_df, 2, var, na.rm = TRUE)
      var_scores[!is.finite(var_scores)] <- 0
      selected_predictors <- names(sort(var_scores, decreasing = TRUE))[1:min(6, ncol(sp_env_df))]
    }
    
    n_predictors_selected <- length(selected_predictors)
    selected_method <- "covsel_filter_only"  # Note: filter only, no embedding
    
    # Apply final selection
    env_selected <- env_masked_subset[[selected_predictors]]
    
    message(sprintf("  Selected %d predictors using %s method", 
                    n_predictors_selected, selected_method))
    
    # Generate background points
    bg_points <- generate_target_group_background(
      env_masked = env_selected,
      accessible_vect = accessible_vect,
      target_species = species_name,
      occ_data = occ_data,
      n_bg_target = 5000
    )
    
    if (is.null(bg_points) || nrow(bg_points) < 10) {
      message("  → Insufficient background points, skipping")
      return(NULL)
    }
    
    # Prepare data for modeling
    model_data <- bind_rows(
      sp_clean %>% mutate(pr_ab = 1),
      bg_points %>% mutate(pr_ab = 0, species = species_name)
    )
    
    # Prepare data for ESM
    coords <- as.matrix(model_data[, c("longitude", "latitude")])
    env_df <- terra::extract(env_selected, coords)
    env_data <- env_df[, -1, drop = FALSE]
    
    # Format data for biomod2/ecospat
    biomod_data <- biomod2::BIOMOD_FormatingData(
      resp.var = model_data$pr_ab,
      expl.var = env_data,
      resp.xy = model_data[, c("longitude", "latitude")],
      resp.name = gsub(" ", "_", species_name)
    )
    
    # Run ESM models
    esm_models <- NULL
    capture.output({
      tryCatch({
        set.seed(123)
        esm_models <- ecospat::ecospat.ESM.Modeling(
          data = biomod_data,
          models = c('GLM', 'GAM', 'MAXNET'),
          NbRunEval = 5,
          DataSplit = 80,
          weighting.score = c("SomersD"),
          parallel = FALSE,
          tune = FALSE
        )
      }, error = function(e) {
        message(sprintf("    ESM failed: %s", e$message))
      })
    })
    
    if (!is.null(esm_models)) {
      # Create ensemble
      thr_D <- if (n_occ < 10) 0.00 else if (n_occ < 15) 0.05 else if (n_occ < 20) 0.10 else 0.15
      
      esm_ensemble <- NULL
      tryCatch({
        esm_ensemble <- ecospat::ecospat.ESM.EnsembleModeling(
          ESM.modeling.output = esm_models,
          weighting.score = c("SomersD"),
          threshold = thr_D
        )
      }, error = function(e) {
        message(sprintf("    ESM ensemble failed: %s", e$message))
      })
      
      # Project to environmental space
      if (!is.null(esm_ensemble)) {
        esm_proj <- NULL
        esm_ensemble_proj <- NULL
        
        tryCatch({
          esm_proj <- ecospat::ecospat.ESM.Projection(
            ESM.modeling.output = esm_models,
            new.env = env_selected
          )
          
          esm_ensemble_proj <- ecospat::ecospat.ESM.EnsembleProjection(
            ESM.prediction.output = esm_proj,
            ESM.EnsembleModeling.output = esm_ensemble
          )
        }, error = function(e) {
          message(sprintf("    ESM projection failed: %s", e$message))
        })
        
        if (!is.null(esm_ensemble_proj)) {
          # Convert to SpatRaster if needed
          if (!inherits(esm_ensemble_proj, "SpatRaster")) {
            esm_ensemble_proj <- terra::rast(esm_ensemble_proj)
          }
          
          # Get layer information
          n_layers <- terra::nlyr(esm_ensemble_proj)
          layer_names <- names(esm_ensemble_proj)
          
          message(sprintf("    ESM projection has %d layers:", n_layers))
          message(sprintf("    Layer names: %s", paste(layer_names, collapse = ", ")))
          
          # Look for EF (Ensemble Forecast) layer - this is the ecospat ensemble!
          if ("EF" %in% layer_names) {
            final_prediction <- esm_ensemble_proj[["EF"]]
            message("    Using EF (Ensemble Forecast) layer")
          } else {
            # Fallback: look for other ensemble patterns
            ensemble_idx <- which(grepl("ensemble|mean|weighted", layer_names, ignore.case = TRUE))
            if (length(ensemble_idx) > 0) {
              final_prediction <- esm_ensemble_proj[[ensemble_idx[length(ensemble_idx)]]]
              message(sprintf("    Using layer: %s", names(final_prediction)))
            } else {
              final_prediction <- esm_ensemble_proj[[n_layers]]
              message(sprintf("    WARNING: Using last layer as fallback: %s", names(final_prediction)))
            }
          }
          
          # Check and normalize values - ROBUST VERSION
          pred_vals <- terra::values(final_prediction, na.rm = TRUE)
          valid_vals <- pred_vals[is.finite(pred_vals)]
          
          if (length(valid_vals) > 0) {
            pred_min <- min(valid_vals, na.rm = TRUE)
            pred_max <- max(valid_vals, na.rm = TRUE)
            
            message(sprintf("    Value range: [%.3f, %.3f] (%d valid cells)", 
                            pred_min, pred_max, length(valid_vals)))
            
            # Normalize if needed
            if (pred_min < 0 || pred_max > 1.1) {
              if (pred_max > pred_min) {
                message("    Normalizing to [0, 1]...")
                final_prediction <- (final_prediction - pred_min) / (pred_max - pred_min)
                
                # Verify
                new_vals <- terra::values(final_prediction, na.rm = TRUE)
                new_valid <- new_vals[is.finite(new_vals)]
                message(sprintf("    After normalization: [%.3f, %.3f]", 
                                min(new_valid, na.rm = TRUE), 
                                max(new_valid, na.rm = TRUE)))
              } else {
                message("    WARNING: All values identical, setting to 0.5")
                final_prediction <- final_prediction * 0 + 0.5
              }
            }
          } else {
            message("    ERROR: No valid prediction values! Skipping species.")
            final_prediction <- NULL
          }
          
          performance_metrics <- list(
            model_type = "ESM",
            n_occurrences = n_occ,
            auc = NA,
            tss = NA,
            n_predictors = n_predictors_selected
          )
        } else {
          final_prediction <- NULL
        }
      } else {
        final_prediction <- NULL
      }
    } else {
      final_prediction <- NULL
    }
    
  } else {
    # Common: Full ensemble modeling with FULL COVSEL (collinearity + embedding)
    message("  Building ensemble model for common species...")
    model_type <- "SDM"
    
    # Predictor selection - FULL COVSEL for common species
    message("  Running full covsel predictor selection (collinearity + embedding)...")
    
    # Generate background points first for covsel
    bg_points_initial <- generate_target_group_background(
      env_masked = env_masked,
      accessible_vect = accessible_vect,
      target_species = species_name,
      occ_data = occ_data,
      n_bg_target = 5000
    )
    
    if (is.null(bg_points_initial) || nrow(bg_points_initial) < 10) {
      message("  → Insufficient background points, skipping")
      return(NULL)
    }
    
    # Prepare pres-bg data for covsel
    pres_bg <- bind_rows(
      sp_clean %>% dplyr::select(longitude, latitude) %>% mutate(pr_ab = 1),
      bg_points_initial %>% dplyr::select(longitude, latitude) %>% mutate(pr_ab = 0)
    )
    
    # Extract environmental data
    xy <- as.matrix(pres_bg[, c("longitude", "latitude")])
    vals <- terra::extract(env_masked, xy)
    vals <- vals[, -1, drop = FALSE]
    covdata_raw <- as.data.frame(lapply(as.data.frame(vals), function(v) as.numeric(v)))
    keep_rows <- stats::complete.cases(covdata_raw) & !is.na(pres_bg$pr_ab)
    covdata_raw <- covdata_raw[keep_rows, , drop = FALSE]
    pa <- pres_bg$pr_ab[keep_rows]
    
    # Remove zero-variance columns
    sd_ok <- vapply(covdata_raw, function(v) {
      s <- stats::sd(v, na.rm = TRUE); is.finite(s) && s > 0
    }, logical(1))
    covdata <- covdata_raw[, sd_ok, drop = FALSE]
    env_work <- env_masked[[which(sd_ok)]]
    
    # Clean names
    colnames(covdata) <- gsub("[^A-Za-z0-9_]", "_", colnames(covdata))
    colnames(covdata) <- make.unique(colnames(covdata), sep = "_")
    names(env_work) <- colnames(covdata)
    
    # Get categories
    cat_vec <- layer_categories[names(env_masked)[which(sd_ok)]]
    cat_vec[is.na(cat_vec)] <- "uncategorized"
    
    # Step A: Collinearity filtering
    message(sprintf("  Starting with %d predictors", ncol(covdata)))
    
    covdata_filt <- tryCatch({
      outA <- covsel::covsel.filter(
        covdata = covdata,
        pa = pa,
        corcut = 0.7,
        categories = cat_vec
      )
      # Handle both possible return types
      if (is.data.frame(outA)) {
        outA
      } else if (is.character(outA)) {
        covdata[, intersect(outA, colnames(covdata)), drop = FALSE]
      } else {
        NULL
      }
    }, error = function(e) {
      message("  covsel.filter failed, using all predictors")
      covdata
    })
    
    if (is.null(covdata_filt) || ncol(covdata_filt) < 2) {
      covdata_filt <- covdata
    }
    message(sprintf("  After Step A (collinearity filtering): %d predictors", ncol(covdata_filt)))
    
    # Step B: Embedded selection
    n_pres <- sum(pa == 1, na.rm = TRUE)
    max_vars <- min(12, max(4, ceiling(log2(max(2, n_pres)))))
    
    covdata_sel <- tryCatch({
      set.seed(42)
      outB <- covsel::covsel.embed(
        covdata = covdata_filt,
        pa = pa,
        algorithms = c("glm","gam","rf"),
        ncov = max_vars,
        maxncov = max_vars,
        nthreads = 1
      )
      # Handle different return types
      if (is.list(outB) && !is.null(outB$covdata)) {
        outB$covdata
      } else if (is.data.frame(outB)) {
        outB
      } else if (is.character(outB)) {
        covdata_filt[, intersect(outB, colnames(covdata_filt)), drop = FALSE]
      } else {
        NULL
      }
    }, error = function(e) {
      message("  covsel.embed failed: ", e$message, " — using top-variance predictors as fallback")
      top_n <- min(max_vars, ncol(covdata_filt))
      vars_by_var <- names(sort(apply(covdata_filt, 2, stats::var, na.rm = TRUE), decreasing = TRUE))[1:top_n]
      covdata_filt[, vars_by_var, drop = FALSE]
    })
    
    if (is.null(covdata_sel) || ncol(covdata_sel) < 2) {
      covdata_sel <- covdata_filt
    }
    
    sel_names <- colnames(covdata_sel)
    message(sprintf("  After Step B (embedding): %d predictors selected", length(sel_names)))
    
    # Map back to raster layers
    matched <- intersect(sel_names, names(env_work))
    if (length(matched) < 2) {
      # Fallback to top variance
      var_scores <- apply(covdata, 2, var, na.rm = TRUE)
      matched <- names(sort(var_scores, decreasing = TRUE))[1:min(6, ncol(covdata))]
    }
    
    env_selected <- env_work[[matched]]
    n_predictors_selected <- terra::nlyr(env_selected)
    selected_method <- "covsel_full"
    selected_predictors <- names(env_selected)
    
    message(sprintf("  Final predictor set: %d variables", n_predictors_selected))
    
    # Generate final background points with selected predictors
    bg_points <- generate_target_group_background(
      env_masked = env_selected,
      accessible_vect = accessible_vect,
      target_species = species_name,
      occ_data = occ_data,
      n_bg_target = 10000
    )
    
    if (is.null(bg_points) || nrow(bg_points) < 10) {
      message("  → Insufficient background points, skipping")
      return(NULL)
    }
    
    # Prepare data for modeling
    model_data <- bind_rows(
      sp_clean %>% mutate(pr_ab = 1),
      bg_points %>% mutate(pr_ab = 0, species = species_name)
    )
    
    # Extract environmental data for modeling
    coords <- as.matrix(model_data[, c("longitude", "latitude")])
    env_df <- terra::extract(env_selected, coords)
    env_data <- env_df[, -1, drop = FALSE]
    
    # Remove NA vs env_selected & keep inside AA
    vals_pb <- terra::extract(env_selected, as.matrix(model_data[, c("longitude", "latitude")]))
    keep_pb <- stats::complete.cases(vals_pb)
    model_data <- model_data[keep_pb, , drop = FALSE]
    env_data <- env_data[keep_pb, , drop = FALSE]
    
    pres_bg_sf_4326 <- sf::st_as_sf(model_data, coords = c("longitude", "latitude"), crs = 4326)
    inside_aa <- sf::st_within(pres_bg_sf_4326, st_as_sf(accessible_vect), sparse = FALSE)[, 1]
    model_data <- model_data[inside_aa, , drop = FALSE]
    env_data <- env_data[inside_aa, , drop = FALSE]
    
    # Deduplicate raster cells
    cell_id <- terra::cellFromXY(env_selected[[1]], as.matrix(model_data[, c("longitude", "latitude")]))
    keep_unique <- !duplicated(cell_id)
    model_data <- model_data[keep_unique, , drop = FALSE]
    env_data <- env_data[keep_unique, , drop = FALSE]
    
    # SPATIAL CV 
    # Convert to sf object in Albers for spatial CV
    model_data_sf <- sf::st_as_sf(model_data, coords = c("longitude", "latitude"), crs = 4326) %>%
      sf::st_transform(ca_albers)
    
    # Project env_selected to Albers for spatial CV
    env_selected_alb <- terra::project(env_selected, ca_albers)
    aa_alb <- sf::st_transform(st_as_sf(accessible_vect), ca_albers)
    env_selected_alb <- terra::mask(terra::crop(env_selected_alb, terra::vect(aa_alb)), terra::vect(aa_alb))
    env_selected_alb <- terra::extend(env_selected_alb, 1)
    
    # Calculate spatial autocorrelation
    autocorr_obj <- blockCV::cv_spatial_autocor(
      x = model_data_sf, 
      column = "pr_ab", 
      r = env_selected_alb, 
      num_sample = 3000, 
      progress = FALSE
    )
    
    # Calculate block size based on autocorrelation
    ext_alb   <- terra::ext(env_selected_alb)
    aa_width  <- ext_alb$xmax - ext_alb$xmin
    aa_height <- ext_alb$ymax - ext_alb$ymin
    cellres   <- mean(terra::res(env_selected_alb))
    raw_range <- suppressWarnings(max(autocorr_obj$range, na.rm = TRUE))
    if (!is.finite(raw_range)) raw_range <- 5 * cellres
    block_size <- max(cellres * 5, min(0.25 * max(aa_width, aa_height), raw_range * 1.25))
    
    # Create spatial blocks
    spatial_blocks <- NULL
    for (k_try in c(5, 4, 3)) {
      attempt <- try(blockCV::cv_spatial(
        x = model_data_sf, 
        column = "pr_ab", 
        r = env_selected_alb,
        k = k_try, 
        size = block_size, 
        selection = "systematic",
        iteration = 50, 
        biomod2 = TRUE, 
        progress = FALSE
      ), silent = TRUE)
      if (!inherits(attempt, "try-error")) { 
        spatial_blocks <- attempt
        break 
      }
    }
    
    cv_strategy <- "spatial_blocks"
    if (is.null(spatial_blocks)) {
      # Fallback to random k-fold if spatial blocking fails
      spatial_blocks <- blockCV::cv_kfold(
        x = model_data_sf, 
        column = "pr_ab", 
        k = 5, 
        biomod2 = TRUE, 
        seed = 42
      )
      cv_strategy <- "random_kfold"
    }
    
    # Create CV table for BIOMOD
    fold_ids <- spatial_blocks$folds_ids
    k <- length(unique(fold_ids))
    cv_tab <- matrix(TRUE, nrow = nrow(model_data), ncol = k)
    for (j in seq_len(k)) {
      cv_tab[, j] <- fold_ids != j
    }
    colnames(cv_tab) <- paste0("_allData_RUN", seq_len(k))
    cv_tab <- cbind(cv_tab, `_allData_allRun` = TRUE)
    storage.mode(cv_tab) <- "logical"
    
    # Prepare BIOMOD data
    biomod_data <- biomod2::BIOMOD_FormatingData(
      resp.var = model_data$pr_ab,
      expl.var = env_data,
      resp.xy = model_data[, c("longitude", "latitude")],
      resp.name = gsub(" ", "_", species_name),
      PA.strategy = NULL
    )
    
    # Create user.val structure for model options
    p <- ncol(env_data)
    user_val <- list(
      RF.binary.randomForest.randomForest = list(
        for_all_datasets = list(do.classif = TRUE, ntree = 800, mtry = max(1, floor(sqrt(p))), nodesize = 5)
      ),
      GBM.binary.gbm.gbm = list(
        for_all_datasets = list(distribution = "bernoulli", n.trees = 800, interaction.depth = 3, 
                                shrinkage = 0.1, bag.fraction = 0.75)
      ),
      MAXNET.binary.maxnet.maxnet = list(
        for_all_datasets = list(classes = "lq", regmult = 1.0)
      ),
      GAM.binary.mgcv.bam = list(
        for_all_datasets = list(select = TRUE, gamma = 1.4)
      )
    )
    
    # Model options using user.defined strategy
    model_opt <- bm_ModelingOptions(
      data.type  = "binary",
      models     = c("GLM", "GAM", "RF", "MAXNET"),
      strategy   = "user.defined",
      user.base  = "bigboss",
      user.val   = user_val,
      bm.format  = biomod_data,
      calib.lines= cv_tab
    )
    
    # Run models with spatial CV
    biomod_models <- biomod2::BIOMOD_Modeling(
      bm.format = biomod_data,
      modeling.id = paste0("CA_", gsub(" ", "_", species_name)),
      models = c("GLM", "GAM", "RF", "MAXNET"),
      OPT.user = model_opt,
      CV.strategy = "user.defined",
      CV.user.table = cv_tab,
      var.import = 2,
      metric.eval = c("TSS", "ROC"),
      CV.do.full.models = TRUE,
      nb.cpu = 1
    )
    
    # Get evaluations and select good models
    eval_df <- as.data.frame(get_evaluations(biomod_models))
    eval_df <- tibble::as_tibble(eval_df)
    
    # Select good models based on TSS
    good_models <- eval_df %>%
      filter(metric.eval == "TSS",
             grepl("_RUN", full.name),
             is.finite(validation),
             validation >= 0.7) %>%
      pull(full.name) %>% unique()
    
    if (length(good_models) < 2) {
      top_algos <- eval_df %>%
        filter(metric.eval == "TSS", grepl("_RUN", full.name)) %>%
        group_by(algo) %>%
        summarize(mean_tss = mean(validation, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_tss)) %>% 
        slice_head(n = 3) %>% 
        pull(algo)
      
      good_models <- eval_df %>%
        filter(metric.eval == "TSS", grepl("_RUN", full.name), algo %in% top_algos) %>%
        pull(full.name) %>% unique()
    }
    
    if (length(good_models) < 1) {
      warning("No models passed TSS threshold; using all available models for ", species_name)
      good_models <- get_built_models(biomod_models)
    }
    
    # Ensemble
    biomod_ensemble <- biomod2::BIOMOD_EnsembleModeling(
      bm.mod = biomod_models,
      models.chosen = good_models,
      em.by = "all",
      em.algo = c("EMmean", "EMwmean"),
      metric.eval = c("TSS", "ROC"),
      var.import = 2,
      EMwmean.decay = 1.6,
      nb.cpu = 1
    )
    
    # Project - FIRST do individual model projection
    sel_names <- biomod_models@expl.var.names
    
    # Get the exact layers that BIOMOD2 used
    # These might be fewer than env_selected if BIOMOD2 dropped some
    ca_env_selected <- env_selected[[sel_names]]  # Select ONLY the layers BIOMOD2 actually used
    ca_env_selected <- terra::mask(terra::crop(ca_env_selected, accessible_vect), accessible_vect)
    
    # The names should already match, but ensure they're correct
    names(ca_env_selected) <- sel_names
    
    current_proj <- BIOMOD_Projection(
      bm.mod = biomod_models,
      proj.name = "current",
      new.env = ca_env_selected,
      models.chosen = good_models,
      metric.binary = "TSS",
      build.clamping.mask = TRUE,
      compress = TRUE,
      nb.cpu = 1
    )
    
    # THEN do ensemble projection using BIOMOD_EnsembleForecasting (NOT BIOMOD_EnsembleProjection)
    ensemble_proj <- BIOMOD_EnsembleForecasting(
      bm.em = biomod_ensemble,
      bm.proj = current_proj,
      models.chosen = "all",
      metric.binary = "TSS",
      compress = TRUE,
      nb.cpu = 1
    )
    
    # Extract predictions
    ensemble_predictions <- get_predictions(ensemble_proj)
    if ("EMwmean" %in% names(ensemble_predictions)) {
      final_prediction <- ensemble_predictions[["EMwmean"]][[1]]
    } else {
      nm_first <- names(ensemble_predictions)[1]
      final_prediction <- ensemble_predictions[[nm_first]][[1]]
    }
    
    # NORMALIZE BIOMOD2 OUTPUT 
    pred_vals <- terra::values(final_prediction, na.rm = TRUE)
    valid_vals <- pred_vals[is.finite(pred_vals)]
    
    if (length(valid_vals) > 0) {
      pred_min <- min(valid_vals, na.rm = TRUE)
      pred_max <- max(valid_vals, na.rm = TRUE)
      
      message(sprintf("    BIOMOD2 value range: [%.3f, %.3f] (%d valid cells)", 
                      pred_min, pred_max, length(valid_vals)))
      
      # Normalize if not already 0-1
      if (pred_min < 0 || pred_max > 1.1) {
        if (pred_max > pred_min) {
          message("    Normalizing BIOMOD2 output to [0, 1]...")
          final_prediction <- (final_prediction - pred_min) / (pred_max - pred_min)
          
          # Verify
          new_vals <- terra::values(final_prediction, na.rm = TRUE)
          new_valid <- new_vals[is.finite(new_vals)]
          message(sprintf("    After normalization: [%.3f, %.3f]", 
                          min(new_valid, na.rm = TRUE), 
                          max(new_valid, na.rm = TRUE)))
        } else {
          message("    WARNING: All BIOMOD2 values identical, setting to 0.5")
          final_prediction <- final_prediction * 0 + 0.5
        }
      } else {
        message("    BIOMOD2 output already in [0, 1] range")
      }
    } else {
      message("    ERROR: No valid BIOMOD2 prediction values! Skipping species.")
      final_prediction <- NULL
    }
    
    # Get performance metrics
    eval_train <- eval_df %>%
      filter(metric.eval == "ROC", full.name == paste0(gsub(" ", "_", species_name), "_allData_allRun")) %>%
      pull(validation) %>% mean(na.rm = TRUE)
    
    eval_test <- eval_df %>%
      filter(metric.eval == "ROC", grepl("_RUN", full.name)) %>%
      pull(validation) %>% mean(na.rm = TRUE)
    
    auc_diff <- eval_train - eval_test
    
    performance_metrics <- list(
      model_type = "SDM",
      n_occurrences = n_occ,
      auc = eval_test,
      tss = mean(eval_df[eval_df$metric.eval == "TSS" & grepl("_RUN", eval_df$full.name), "validation"], na.rm = TRUE),
      n_predictors = n_predictors_selected,
      cv_strategy = cv_strategy,
      auc_diff = auc_diff
    )
  }
  
  # Step 9: Apply MESS filtering (skip for ultra-rare) - FIXED
  if (n_occ >= 7 && !is.null(final_prediction)) {
    message("  Applying MESS filtering to prevent extrapolation...")
    
    # Extract training environmental values (presences only)
    training_coords <- model_data[model_data$pr_ab == 1, c("longitude", "latitude")]
    training_env <- terra::extract(env_selected, training_coords)
    training_env <- training_env[, -1, drop = FALSE]
    
    # Convert SpatRaster to RasterStack for dismo::mess
    env_raster <- raster::stack(env_selected)
    
    # Calculate MESS
    mess_map <- tryCatch({
      dismo::mess(env_raster, training_env, full = TRUE)
    }, error = function(e) {
      message("    MESS calculation failed: ", e$message)
      NULL
    })
    
    if (!is.null(mess_map)) {
      # Convert back to SpatRaster
      mess_spatraster <- terra::rast(mess_map)
      
      # CRITICAL: Extract ONLY the MESS layer (usually first layer)
      # dismo::mess with full=TRUE returns: layer 1 = MESS, layers 2+ = MoD for each variable
      mess_layer <- mess_spatraster[[1]]  # ONLY the MESS index
      
      # Create binary mask (MESS >= 0 means interpolation)
      mess_binary <- mess_layer >= 0
      
      # Apply mask to prediction (NOT replace final_prediction!)
      final_prediction <- final_prediction * mess_binary
      
      message("    MESS filtering applied: removed extrapolation areas")
    }
  }

  # Step 10: Apply dispersal constraints (OPTIMIZED & ROBUST)
  dispersal_constraint_applied <- FALSE
  dispersal_constraint_metadata <- NULL
  
  if (n_occ >= 7 && n_occ <= 100 && !is.null(final_prediction)) {
    
    coords_matrix <- as.matrix(sp_clean[, c("longitude", "latitude")])
    
    if (nrow(coords_matrix) > 1) {
      sp_pts_wgs <- sp::SpatialPoints(coords_matrix, proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
      sp_pts_alb <- sp::spTransform(sp_pts_wgs, sp::CRS(ca_albers))
      
      dist_matrix <- sp::spDists(sp_pts_alb)
      diag(dist_matrix) <- NA
      nn_distances <- apply(dist_matrix, 1, min, na.rm = TRUE)
      max_nn_dist <- max(nn_distances, na.rm = TRUE)
      median_nn_dist <- median(nn_distances, na.rm = TRUE)
      clustering_ratio <- max_nn_dist / median_nn_dist
      
      if (clustering_ratio > 5) {
        message(sprintf("  Applying dispersal constraints (clustering ratio=%.2f)...", clustering_ratio))
        
        if (terra::nlyr(final_prediction) > 1) final_prediction <- final_prediction[[1]]
        
        sp_clean_with_pr <- if (!"pr_ab" %in% names(sp_clean)) {
          mutate(sp_clean, pr_ab = 1)
        } else sp_clean
        
        occ_preds <- terra::extract(final_prediction, sp_clean_with_pr[, c("longitude", "latitude")])[, 2]
        
        if (exists("bg_points") && !is.null(bg_points) && nrow(bg_points) > 0) {
          bg_preds <- terra::extract(final_prediction, bg_points[, c("longitude", "latitude")])[, 2]
          thr_ms <- calculate_OR10_threshold(
            obs = c(rep(1, length(occ_preds)), rep(0, length(bg_preds))),
            pred = c(occ_preds, bg_preds)
          )
        } else {
          thr_ms <- median(occ_preds, na.rm = TRUE)
        }
        
        thr_ms <- max(as.numeric(thr_ms), 0.01, na.rm = TRUE)
        cont_rlay <- raster::raster(final_prediction)
        
        # ROBUST COERCION HELPER
        as_spatraster <- function(x) {
          if (is.null(x)) return(NULL)
          if (inherits(x, "SpatRaster")) return(x)
          if (inherits(x, "Raster")) return(terra::rast(x))
          if (is.list(x)) {
            cand <- x$cont_suit %||% x$cont.suit %||% x$suitability %||% x$cont
            if (!is.null(cand)) return(as_spatraster(cand))
          }
          NULL
        }
        
        constrained_result <- tryCatch({
          suppressWarnings(
            flexsdm::msdm_posteriori(
              records = sp_clean_with_pr,
              x = "longitude",
              y = "latitude",
              pr_ab = "pr_ab",
              method = "obr",
              cont_suit = cont_rlay,
              thr = thr_ms
            )
          )
        }, error = function(e) {
          message(sprintf("    Dispersal constraint failed: %s", e$message))
          NULL
        })
        
        constrained_sr <- as_spatraster(constrained_result)
        
        if (!is.null(constrained_sr)) {
          final_prediction <- constrained_sr
          dispersal_constraint_applied <- TRUE
          message("    Dispersal constraints applied successfully")
        }
        
        dispersal_constraint_metadata <- list(
          method = ifelse(dispersal_constraint_applied, "obr", "failed"),
          clustering_ratio = clustering_ratio,
          max_nn_distance_m = max_nn_dist
        )
        
      } else {
        message(sprintf("  Skipping dispersal constraints: not sufficiently clustered (ratio=%.2f)", clustering_ratio))
        dispersal_constraint_metadata <- list(method = "skipped_not_clustered", clustering_ratio = clustering_ratio)
      }
    }
  } else if (n_occ >= 7) {
    message("  Skipping dispersal constraints: widespread species (n>100)")
  }
  
  
  # Step 11: Apply strict ecoregion masking (no buffer)
  if (!is.null(final_prediction)) {
    message("  Applying strict ecoregion masking...")
    final_prediction <- terra::mask(
      terra::crop(final_prediction, strict_ecoregion_vect),
      strict_ecoregion_vect
    )
  }
  
  # Step 12: Align to template
  if (!is.null(final_prediction)) {
    final_prediction_aligned <- terra::resample(final_prediction, template, method = "bilinear")
    final_prediction_aligned <- terra::mask(final_prediction_aligned, template)
    # Ensure single layer with simple name
    if (terra::nlyr(final_prediction_aligned) > 1) {
      final_prediction_aligned <- final_prediction_aligned[[1]]
    }
    names(final_prediction_aligned) <- species_name
  } else {
    final_prediction_aligned <- NULL
  }
  
  # Step 13: Calculate threshold for binary map
  binary_prediction <- NULL
  threshold_value <- NA
  
  if (!is.null(final_prediction_aligned)) {
    if (n_occ < 7) {
      # For ultra-rare species, create binary from buffer/hull
      # Option 1: Simple presence within hull (anywhere > 0)
      threshold_value <- 0.01  # or use 0 if you want strict hull boundary
      binary_prediction <- final_prediction_aligned > threshold_value
      
      # Option 2 (alternative): Use median of the distance decay as threshold
      # threshold_value <- terra::global(final_prediction_aligned, "median", na.rm = TRUE)[1,1]
      # binary_prediction <- final_prediction_aligned >= threshold_value
      
      names(binary_prediction) <- species_name
      
    } else {
      # Existing code for n >= 7 species
      # Handle bg_points for each branch
      if (!exists("bg_points")) {
        # For ultra-rare, no background points
        bg_points <- data.frame(longitude = numeric(0), latitude = numeric(0))
      }
      
      # Extract predictions at occurrence points
      occ_preds <- terra::extract(final_prediction_aligned, sp_clean[, c("longitude", "latitude")])[, 2]
      
      if (nrow(bg_points) > 0) {
        bg_preds <- terra::extract(final_prediction_aligned, bg_points[, c("longitude", "latitude")])[, 2]
        
        # Calculate threshold
        threshold_value <- calculate_OR10_threshold(
          obs = c(rep(1, length(occ_preds)), rep(0, length(bg_preds))),
          pred = c(occ_preds, bg_preds)
        )
      } else {
        # No background points, use median of occurrence predictions
        threshold_value <- median(occ_preds, na.rm = TRUE)
      }
      
      # Create binary prediction with simple name
      binary_prediction <- final_prediction_aligned >= threshold_value
      names(binary_prediction) <- species_name
    }
  }
  
  # === Step 14: Write each SDM as an individual file (no stacking) ===
  if (!is.null(final_prediction_aligned)) {
    taxon_dir <- file.path(taxon_output_base, gsub("[^[:alnum:]_-]", "_", species_taxon))
    cont_dir  <- file.path(taxon_dir, "continuous")
    dir.create(cont_dir, showWarnings = FALSE, recursive = TRUE)
    
    cont_path <- file.path(cont_dir, paste0(gsub("[^[:alnum:]_-]", "_", species_name), "_continuous.tif"))
    terra::writeRaster(final_prediction_aligned, cont_path, overwrite = TRUE)
    
    # Append manifest row
    cont_manifest_path <- file.path(cont_dir, paste0(species_taxon, "_continuous_manifest.csv"))
    cont_row <- tibble(
      species = species_name,
      taxon = species_taxon,
      model_type = model_type,
      n_occurrences_raw = n_occ_raw,
      n_occurrences_final = n_occ,
      n_predictors = n_predictors_selected,
      predictor_selection_method = selected_method,
      selected_predictors = paste(selected_predictors, collapse = ", "),
      threshold_method = ifelse(!is.na(threshold_value), THRESHOLD_METHOD, NA),
      threshold_value = threshold_value,
      auc = performance_metrics$auc %||% NA,
      tss = performance_metrics$tss %||% NA,
      ecoregions = paste(unique_ecoregions, collapse = ", ")
    )
    write.table(cont_row, cont_manifest_path,
                sep = ",", row.names = FALSE, col.names = !file.exists(cont_manifest_path),
                append = TRUE)
    
    rm(final_prediction_aligned)
    gc()
  }
  
  if (!is.null(binary_prediction)) {
    taxon_dir <- file.path(taxon_output_base, gsub("[^[:alnum:]_-]", "_", species_taxon))
    bin_dir   <- file.path(taxon_dir, "binary")
    dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)
    
    bin_path <- file.path(bin_dir, paste0(gsub("[^[:alnum:]_-]", "_", species_name), "_binary.tif"))
    terra::writeRaster(binary_prediction, bin_path, overwrite = TRUE)
    
    # Append manifest row
    bin_manifest_path <- file.path(bin_dir, paste0(species_taxon, "_binary_manifest.csv"))
    bin_row <- tibble(
      species = species_name,
      taxon = species_taxon,
      model_type = model_type,
      n_occurrences_raw = n_occ_raw,
      n_occurrences_final = n_occ,
      n_predictors = n_predictors_selected,
      predictor_selection_method = selected_method,
      selected_predictors = paste(selected_predictors, collapse = ", "),
      threshold_method = ifelse(!is.na(threshold_value), THRESHOLD_METHOD, NA),
      threshold_value = threshold_value,
      auc = performance_metrics$auc %||% NA,
      tss = performance_metrics$tss %||% NA,
      ecoregions = paste(unique_ecoregions, collapse = ", ")
    )
    write.table(bin_row, bin_manifest_path,
                sep = ",", row.names = FALSE, col.names = !file.exists(bin_manifest_path),
                append = TRUE)
    
    rm(binary_prediction)
    gc()
  }
  
}
  
#### Main processing loop ####

# Base species (2-word names)
base_species <- unique(occ_data$species)

# Subspecies (3-word names only)
subspecies <- occ_data %>%
  filter(str_count(species_full, "\\S+") == 3) %>%  # keep only 3-word names
  pull(species_full) %>%
  unique()

# # To run species and subspecies
unique_species <- unique(c(base_species, subspecies))

#To only run subspecies
# unique_species <- unique(subspecies)

message(sprintf("Found %d total modeling targets (%d species, %d subspecies)",
                length(unique_species), length(base_species), length(subspecies)))

message(sprintf("\n=== Starting SDM processing for %d species ===\n", length(unique_species)))

## Sequential run — no accumulation in memory + timer
run_silent <- function(expr) {
  suppressMessages(suppressWarnings(capture.output(
    invisible(force(expr))
  )))
}

n_species <- length(unique_species)
failed_species <- character(0)  # Track failures

# Create progress bar with elapsed/ETA fields
pb <- progress_bar$new(
  format = "[:bar] :percent | :current/:total species | Elapsed: :elapsed | ETA: :eta",
  total = n_species,
  width = 70,
  clear = FALSE
)

start_time <- Sys.time()

for (i in seq_along(unique_species)) {
  sp <- unique_species[i]
  
  tryCatch({
    run_silent({
      sdm_for_species(
        species_name = sp,
        occ_data = occ_data,
        env_stack = env_stack,
        layer_categories = layer_categories,
        category_reference = category_reference,
        jepson_proj = jepson_proj,
        ca_boundary_proj = ca_boundary_proj,
        biomod_output_dir = biomod_output_dir
      )
    })
  }, error = function(e) {
    failed_species <<- c(failed_species, sp)
    message(sprintf("✗ Failed: %s", sp))
  })
  
  pb$tick()
  gc()
}

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat(sprintf("\n✅ SDM processing complete for %d species in %.1f minutes (%.2f hours)\n",
            n_species, as.numeric(total_time), as.numeric(total_time) / 60))
cat(sprintf("   Successful: %d | Failed: %d\n", n_species - length(failed_species), length(failed_species)))

# Write failed species log
if (length(failed_species) > 0) {
  failed_log <- file.path(output_dir, "failed_species.csv")
  write.csv(data.frame(species = failed_species), failed_log, row.names = FALSE)
  cat(sprintf("   Failed species logged to: %s\n", failed_log))
}

#### Copy SDM outputs to Data_Calscape ####
## NOTE - CHRIS can you add a few lines here that will generate the final SDMs_Calscape folder we send to Calscape here?
# It should be saved at "Data_Calscape/SDMs/SDMs_Calscape"
# Please also join and make a copy of sdm_data.csv at "Data_Calscape/SDMs/SDMs_Inputs", which should be currently in "Data_Clean/SDMs/SDM_Inputs/occ_split"


#### Plotting and Troubleshooting ####

# #### Plot all SDMs with occurrence points (combined across taxa) ####
# 
# # --- Auto California outline ---
# ca_outline <- ne_states(country = "United States of America", returnclass = "sf") |>
#   filter(name == "California") |> vect()
# 
# # --- Base directories ---
# taxon_base <- file.path(output_dir, "sdm_by_taxon")
# plot_dir <- file.path(output_dir, "plots_all_taxa")
# dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # --- Helper to get occurrence coordinates per species ---
# get_occ_points <- function(sp, crs_target) {
#   # Match either straight species (2-word) or subspecies (3-word)
#   sub <- occ_data[occ_data$species == sp | occ_data$species_full == sp,
#                   c("longitude", "latitude")]
# 
#   if (nrow(sub) == 0) return(NULL)
# 
#   pts <- st_as_sf(sub, coords = c("longitude", "latitude"), crs = 4326)
#   pts <- st_transform(pts, crs_target)
#   st_coordinates(pts)
# }
# 
# # --- Quick plotting function ---
# plot_sdm <- function(r_layer, bin_layer, sp_name, ca_outline, occ_xy, out_dir) {
#   sp_clean <- gsub("[^[:alnum:]_\\-]+", "_", sp_name)
# 
#   # Continuous
#   png(file.path(out_dir, paste0(sp_clean, "_continuous.png")), width = 1200, height = 1000, res = 150)
#   plot(r_layer, col = terrain.colors(100)[100:1], main = sp_name, legend = TRUE, axes = FALSE, box = FALSE)
#   plot(ca_outline, add = TRUE, border = "black", lwd = 1.5)
#   if (!is.null(occ_xy)) points(occ_xy, pch = 19, col = "black", cex = 0.5)
#   dev.off()
# 
#   # Binary
#   png(file.path(out_dir, paste0(sp_clean, "_binary.png")), width = 1200, height = 1000, res = 150)
#   plot(bin_layer, col = c("white", "lightgreen"), main = paste(sp_name, "(binary)"),
#        legend = FALSE, axes = FALSE, box = FALSE)
#   plot(ca_outline, add = TRUE, border = "black", lwd = 1.5)
#   if (!is.null(occ_xy)) points(occ_xy, pch = 19, col = "black", cex = 0.5)
#   legend("bottomleft", legend = c("Unsuitable", "Suitable"), fill = c("white", "lightgreen"), bty = "n")
#   dev.off()
# }
# 
# # --- Gather all continuous and binary SDMs across taxa ---
# cont_files <- list.files(taxon_base, pattern = "_continuous\\.tif$", recursive = TRUE, full.names = TRUE)
# bin_files  <- list.files(taxon_base, pattern = "_binary\\.tif$", recursive = TRUE, full.names = TRUE)
# 
# message(sprintf("\nFound %d continuous and %d binary SDM files to plot", length(cont_files), length(bin_files)))
# 
# # --- Iterate over all continuous files (individual species files) ---
# for (cf in cont_files) {
#   # Extract species name from filename
#   # Example: "Abies_concolor_continuous.tif" -> "Abies concolor"
#   sp_clean_from_file <- tools::file_path_sans_ext(basename(cf))
#   sp_name <- gsub("_continuous$", "", sp_clean_from_file)
#   sp_name <- gsub("_", " ", sp_name)  # Convert underscores back to spaces
# 
#   taxon_name <- basename(dirname(dirname(cf)))  # parent taxon folder
# 
#   # Load continuous raster
#   r_cont <- terra::rast(cf)
# 
#   # Find corresponding binary file
#   bf <- gsub("/continuous/", "/binary/", cf)
#   bf <- gsub("_continuous\\.tif$", "_binary.tif", bf)
#   r_bin <- if (file.exists(bf)) {
#     terra::rast(bf)
#   } else {
#     # fallback threshold of 0.5 on continuous raster
#     message(sprintf("  Warning: No binary file for %s, using 0.5 threshold", sp_name))
#     r_cont >= 0.5
#   }
# 
#   # Get occurrence points
#   occ_xy <- get_occ_points(sp_name, terra::crs(r_cont))
# 
#   # Plot
#   plot_sdm(r_cont, r_bin, sp_name, ca_outline, occ_xy, plot_dir)
#   message(sprintf("  Saved plots for %s (taxon: %s)", sp_name, taxon_name))
# }
# 
# message("\n=== All SDM plots (continuous + binary) saved to: ===")
# message(plot_dir)
# 
# #### Check missing species ####
# 
# # 1) Ground truth: all target species
# target_species <- sort(unique(occ_data$species))
# 
# # 2) Collect species names from files written to disk
# taxon_base <- file.path(output_dir, "sdm_by_taxon")
# 
# cont_files <- list.files(taxon_base, pattern = "_continuous\\.tif$",
#                          recursive = TRUE, full.names = TRUE)
# bin_files  <- list.files(taxon_base, pattern = "_binary\\.tif$",
#                          recursive = TRUE, full.names = TRUE)
# 
# # Helper to extract species names from filenames
# get_species_names_from_files <- function(files, pattern) {
#   if (!length(files)) return(character(0))
#   basenames <- basename(files)
#   # Remove the pattern suffix and convert underscores to spaces
#   species <- gsub(pattern, "", basenames)
#   species <- gsub("_", " ", species)
#   unique(species)
# }
# 
# cont_species <- sort(get_species_names_from_files(cont_files, "_continuous\\.tif$"))
# bin_species  <- sort(get_species_names_from_files(bin_files, "_binary\\.tif$"))
# 
# # 3) Who's missing?
# missing_in_cont <- setdiff(target_species, cont_species)
# missing_in_bin  <- setdiff(target_species, bin_species)
# missing_in_either <- setdiff(target_species, union(cont_species, bin_species))
# 
# # 4) Quick summaries
# cat("\n=== Coverage summary ===\n")
# cat(sprintf("Targets: %d\n", length(target_species)))
# cat(sprintf("In continuous: %d\n", length(intersect(target_species, cont_species))))
# cat(sprintf("In binary:     %d\n", length(intersect(target_species, bin_species))))
# cat(sprintf("Missing (continuous): %d\n", length(missing_in_cont)))
# cat(sprintf("Missing (binary):     %d\n", length(missing_in_bin)))
# cat(sprintf("Missing (both):       %d\n", length(missing_in_either)), "\n")
# 
# # 5) Optional: write lists
# write.csv(data.frame(species = missing_in_cont),
#           file.path(output_dir, "missing_in_continuous.csv"), row.names = FALSE)
# write.csv(data.frame(species = missing_in_bin),
#           file.path(output_dir, "missing_in_binary.csv"), row.names = FALSE)
# write.csv(data.frame(species = missing_in_either),
#           file.path(output_dir, "missing_in_either.csv"), row.names = FALSE)
