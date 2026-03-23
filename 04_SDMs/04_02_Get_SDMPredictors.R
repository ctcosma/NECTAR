# =========================
# California SDM predictor stack (EPSG:4326, 5 arcmin)
# End-to-end: fetch, derive, align, and write
# =========================

# Packages
required <- c(
  "terra","sf","dplyr","geodata",
  "rnaturalearth","rnaturalearthdata",
  "nhdplusTools","FedData"
)
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(terra)
library(sf)
library(dplyr)
library(geodata)
library(rnaturalearth)
library(rnaturalearthdata)
library(nhdplusTools)
library(FedData)

# ---- Paths & options ----
data_dir <- file.path(getwd(), "Data_Clean/SDMs/predictors_cache")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
geodata_path(data_dir)  # cache for geodata downloads

out_dir <- file.path(getwd(), "Data_Clean/SDMs/predictors_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Area of Interest (California + 100km buffer) ----
ca_sf <- rnaturalearth::ne_states(country = "United States of America",
                                  returnclass = "sf") |>
  filter(name == "California") |>
  st_make_valid() |>
  st_transform("EPSG:4326")

# Create 100km buffer in distance-preserving projection
ca_sf <- st_transform(ca_sf, "EPSG:3310")  # California Albers Equal Area
ca_sf <- st_buffer(ca_sf, dist = 100000)  # 100km = 100,000 meters

# Transform buffer back to EPSG:4326
ca_sf <- st_transform(ca_sf, "EPSG:4326")

# Convert both to SpatVector for terra operations
ca_v <- vect(ca_sf)

#outline of just ca to compare
ca_outline <- rnaturalearth::ne_states(country = "United States of America",
                                  returnclass = "sf") |>
  filter(name == "California") |>
  st_make_valid() |>
  st_transform("EPSG:4326")

# ---- Target grid: EPSG:4326, 5 arcmin ----
# We'll anchor all resampling to the WorldClim BIO grid (native 5 arcmin), then align others to it.
bio5 <- geodata::worldclim_global(var = "bio", res = 5, path = data_dir)
bio5_ca <- crop(bio5, ca_v)
bio5_ca <- mask(bio5_ca, ca_v, touches = TRUE)
template <- bio5_ca[[1]]  # now template exists

# Plot raster
plot(template)
# Overlay California outline in red
plot(st_geometry(ca_outline), add = TRUE, border = "red", lwd = 2)

# # ---- Climate (PRISM normals) ----
# # Using 30-year (1991–2020) monthly normals from PRISM (4 km)
# # Reference: PRISM Climate Group, Oregon State University (https://prism.oregonstate.edu)
# 
# # Install and load if not already present
# if (!requireNamespace("prism", quietly = TRUE)) install.packages("prism")
# library(prism)
# 
# # Set PRISM data directory (use same cache path as other predictors)
# options(prism.path = data_dir)
# 
# # Download monthly mean temperature (tmean) and precipitation normals
# prism::get_prism_normals(type = "tmean", resolution = "4km", annual = FALSE, keepZip = FALSE)
# prism::get_prism_normals(type = "ppt",   resolution = "4km", annual = FALSE, keepZip = FALSE)
# 
# # List local PRISM rasters
# tmean_files <- prism::ls_prism_data(name = TRUE) |> grep("tmean", ., value = TRUE)
# ppt_files   <- prism::ls_prism_data(name = TRUE) |> grep("ppt",   ., value = TRUE)
# 
# # Load all months as SpatRaster stacks
# tmean_stack <- terra::rast(prism::pd_to_file(tmean_files))
# ppt_stack   <- terra::rast(prism::pd_to_file(ppt_files))
# 
# # Reproject to EPSG:4326 for consistency
# tmean_stack <- terra::project(tmean_stack, "EPSG:4326", method = "bilinear")
# ppt_stack   <- terra::project(ppt_stack,   "EPSG:4326", method = "bilinear")
# 
# # Crop/mask to California
# tmean_ca <- crop(tmean_stack, ca_v) |> mask(ca_v)
# ppt_ca   <- crop(ppt_stack,   ca_v) |> mask(ca_v)
# 
# # ---- Derive simple bioclimatic variables from PRISM ----
# # These replicate BIO1, BIO12, BIO5, BIO6, BIO7, BIO4-style metrics using monthly normals
# 
# # Annual mean temperature (BIO1)
# bio1 <- mean(tmean_ca, na.rm = TRUE)
# 
# # Mean diurnal range not available directly — skip or compute from tmax/tmin if desired
# # Annual precipitation (BIO12)
# bio12 <- sum(ppt_ca, na.rm = TRUE)
# 
# # Max temperature of warmest month (BIO5 proxy)
# bio5 <- max(tmean_ca, na.rm = TRUE)
# 
# # Min temperature of coldest month (BIO6 proxy)
# bio6 <- min(tmean_ca, na.rm = TRUE)
# 
# # Temperature annual range (BIO7 = BIO5 - BIO6)
# bio7 <- bio5 - bio6
# 
# # Precipitation seasonality (BIO15 proxy; coefficient of variation)
# bio15 <- (sd(ppt_ca, na.rm = TRUE) / mean(ppt_ca, na.rm = TRUE)) * 100
# 
# # Combine to one SpatRaster stack
# bio_prism <- c(bio1, bio12, bio5, bio6, bio7, bio15)
# names(bio_prism) <- c("bio1_tmean_annual",
#                       "bio12_ppt_annual",
#                       "bio5_tmean_max_month",
#                       "bio6_tmean_min_month",
#                       "bio7_tmean_range",
#                       "bio15_ppt_seasonality")
# 
# # ---- Create 5-arcmin template aligned to PRISM extent ----
# # Generate template from PRISM grid (≈4 km native)
# template <- terra::rast(bio1)
# template <- terra::resample(template, bio1, method = "near")  # Keep structure
# template <- project(template, "EPSG:4326")                   # ensure lat/lon
# 
# # Aggregate PRISM (4 km → ~5 arcmin ≈ 9 km) for integration with other predictors
# bio_prism_agg <- terra::aggregate(bio_prism, fact = 2, fun = mean, na.rm = TRUE)
# bio_prism_agg <- resample(bio_prism_agg, template, method = "bilinear")
# 
# # Final crop/mask for California extent
# bio_prism_ca <- crop(bio_prism_agg, ca_v) |> mask(ca_v)

# ---- Topography (SRTM 30") -> derived terrain, aggregated to 5' ----
# Elevation at ~1 km (30") then derive slope/aspect/roughness and aggregate to template
elev_us <- geodata::elevation_30s(country = "United States of America", path = data_dir)
elev_mx <- geodata::elevation_30s(country = "Mexico", path = data_dir)

# Merge the two rasters
elev_30s <- merge(elev_us, elev_mx)
elev_ca  <- crop(elev_30s, ca_v) |> mask(ca_v, touches = TRUE)
# Derive terrain on native resolution for better geomorphometry, then aggregate to 5'
slope_30s <- terrain(elev_ca, v = "slope", unit = "degrees")
aspect_30s<- terrain(elev_ca, v = "aspect", unit = "degrees")
tri_30s   <- terrain(elev_ca, v = "TRI")

# Aggregate to ~5 arcmin (~10x10 of 30" cells ≈ 10x10=100; use mean for continuous, circular for aspect)
agg_factor_x <- round((res(template)[1]) / (res(elev_ca)[1]))   # ~10
agg_factor_y <- round((res(template)[2]) / (res(elev_ca)[2]))   # ~10
agg_mean <- function(x, ...){ mean(x, na.rm = TRUE) }

elev_5m  <- aggregate(elev_ca, fact = c(agg_factor_x, agg_factor_y), fun = agg_mean)
slope_5m <- aggregate(slope_30s, fact = c(agg_factor_x, agg_factor_y), fun = agg_mean)
tri_5m   <- aggregate(tri_30s,   fact = c(agg_factor_x, agg_factor_y), fun = agg_mean)

# Aspect: convert to sin/cos components before aggregation to avoid circularity, then back to degrees
asp_rad  <- (aspect_30s * pi) / 180
asp_sin  <- sin(asp_rad); asp_cos <- cos(asp_rad)
asp_sin5 <- aggregate(asp_sin, fact = c(agg_factor_x, agg_factor_y), fun = agg_mean)
asp_cos5 <- aggregate(asp_cos, fact = c(agg_factor_x, agg_factor_y), fun = agg_mean)
aspect_5m<- atan2(asp_sin5, asp_cos5) * 180 / pi
aspect_5m <- (aspect_5m + 360) %% 360  # 0–360

# Align to WorldClim 5' grid exactly (handles any small rounding)
elev_5m   <- resample(elev_5m,   template, method = "bilinear")
slope_5m  <- resample(slope_5m,  template, method = "bilinear")
aspect_5m <- resample(aspect_5m, template, method = "bilinear")
tri_5m    <- resample(tri_5m,    template, method = "bilinear")

names(elev_5m)   <- "elev"
names(slope_5m)  <- "slope"
names(aspect_5m) <- "aspect"
names(tri_5m)    <- "tri"

# ---- Land cover fractions (global, ~30") and Human Footprint (global) ----
# Fractions of built/trees/cropland/etc. + footprint (anthropogenic pressure)
lc_vars <- c("built","trees","shrubs","grassland","cropland","bare","water","wetland")
lc_list <- lapply(lc_vars, function(v) {
  r <- geodata::landcover(var = v, path = data_dir)           # returns 0–1 fraction
  r <- crop(r, ca_v) 
  r <- mask(r, ca_v, touches = TRUE)
  resample(r, template, method = "bilinear")
})
lc <- do.call(c, lc_list)
names(lc) <- paste0("lc_", lc_vars)

# # ---- Landscape Intactness (aligned to 5' EPSG:4326 template) ----
# 
# landscape_intactness_gdb <- '/Users/chriscosma/Desktop/CBI/Current Projects/2025_CDFA:CNPS/Calscape SDMs/input/environmental/Landscape Intactness (1 km), California - 2025/data/p20/data.gdb'
# 
# # Load polygon layer
# landscape_intactness_vec <- terra::vect(
#   landscape_intactness_gdb,
#   layer = "Landscape_Intactness__1_km___California___2025_Update_wWGRtw"
# )
# 
# # Numeric attribute to rasterize
# value_field <- "High_Intactness"
# 
# # 1. Reproject polygons to EPSG:4326 to match your final template
# landscape_intactness_vec <- terra::project(landscape_intactness_vec, "EPSG:4326")
# 
# # 2. Clean geometry just in case
# landscape_intactness_vec <- terra::makeValid(landscape_intactness_vec)
# 
# # 3. Rasterize directly to the 5-arcmin WorldClim template grid
# #    This guarantees identical extent/resolution/CRS to all other layers.
# landscape_intactness_5m <- terra::rasterize(
#   landscape_intactness_vec,
#   template,                 # EPSG:4326, 5'
#   field = value_field,
#   fun   = "mean",
#   background = NA
# )
# 
# # 4. Crop and mask to California boundary
# landscape_intactness_5m <- crop(landscape_intactness_5m, ca_v)
# landscape_intactness_5m <- mask(landscape_intactness_5m, ca_v, touches = TRUE)
# 
# # 5. Standardize name
# names(landscape_intactness_5m) <- "landscape_intactness"

#### Human footprint ####

# From: https://www.nature.com/articles/s41597-022-01284-8

hfp_path <- 'Data_Clean/SDMs/SDM_Inputs/environmental/human_footprint/hfp2018.tif'

# Load the raster
hfp_raw <- terra::rast(hfp_path)

# Check and reproject if needed
if (!identical(terra::crs(hfp_raw), terra::crs(ca_v))) {
  message("  Reprojecting Human Footprint to EPSG:4326...")
  hfp_raw <- terra::project(hfp_raw, "EPSG:4326", method = "bilinear")
}

# Crop to California + 100km buffer extent
hfp_cropped <- terra::crop(hfp_raw, ca_v)

# Mask to the buffer area
hfp_masked <- terra::mask(hfp_cropped, ca_v, touches = TRUE)

# Resample to match your 5 arcmin template
hfp_5m <- terra::resample(hfp_masked, template, method = "bilinear")

# Name it
names(hfp_5m) <- "human_footprint"

# # Optional: Annual NLCD (CONUS) + Fractional Impervious (if you prefer NLCD coding)
# # Uses "Annual NLCD" endpoint (recommended; 2021 L48 on legacy server is not available)
# # Comment out if not needed.
# nlcd_lc  <- try({
#   FedData::get_nlcd_annual(template = ca_sf, label = "CA", year = 2023,
#                            product = "LndCov", region = "CU")
# }, silent = TRUE)
# nlcd_imp <- try({
#   FedData::get_nlcd_annual(template = ca_sf, label = "CA", year = 2023,
#                            product = "FctImp", region = "CU")
# }, silent = TRUE)
# 
# if (inherits(nlcd_lc, "SpatRaster")) {
#   nlcd_mode_5m <- resample(nlcd_lc, template, method = "near")        # dominant class per 5' cell
#   names(nlcd_mode_5m) <- "nlcd2023_mode"
# } else {
#   nlcd_mode_5m <- NULL
# }
# 
# if (inherits(nlcd_imp, "SpatRaster")) {
#   nlcd_imp_5m <- resample(nlcd_imp/100, template, method = "bilinear") # 0–1 fraction
#   names(nlcd_imp_5m) <- "nlcd2023_frac_impervious"
# } else {
#   nlcd_imp_5m <- NULL
# }

# ---- Soils (SoilGrids-derived, via geodata::soil_world) ----
# Variables at standard depths; compute 0–30 cm weighted means.
# All use 0–5, 5–15, 15–30 cm except "clay", which skips the problematic 0–5 cm layer.
soil_vars <- c("sand","silt","clay","phh2o","soc","bdod")

soil_stack <- list()

for (sv in soil_vars) {
  if (sv == "clay") {
    depths <- c(5, 15, 30)[-1]     # use only 5–15 and 15–30 cm for clay
    thick  <- c(10, 15)            # total 25 cm
  } else {
    depths <- c(5, 15, 30)         # use 0–5, 5–15, and 15–30 cm
    thick  <- c(5, 10, 15)         # total 30 cm
  }
  
  layers <- lapply(depths, function(dpt) {
    r <- geodata::soil_world(var = sv, depth = dpt, path = data_dir)
    crop(r, ca_v) |> mask(ca_v, touches = TRUE)
  })
  
  # Weighted mean by layer thickness
  num <- 0; den <- 0
  for (i in seq_along(layers)) {
    li <- layers[[i]]
    li <- resample(li, layers[[1]], method = "bilinear")
    num <- num + (li * thick[i]); den <- den + thick[i]
  }
  
  wmean <- num / den
  soil_stack[[sv]] <- resample(wmean, template, method = "bilinear")
}

soil <- rast(soil_stack)
names(soil) <- paste0("soil_", soil_vars, "_0_30cm")

# ---- Population density (GPW, 5 arcmin) ----
pop5 <- geodata::population(year = 2020, res = 5, path = data_dir)
pop5_ca <- crop(pop5, ca_v)
pop5_ca <- mask(pop5_ca, ca_v, touches = TRUE)
pop_5m  <- resample(pop5_ca, template, method = "bilinear")
names(pop_5m) <- "pop2020_gpw"

# ---- Distance to coast (m) and distance to streams (m), projected compute then back to EPSG:4326 ----
# Compute distances in California Albers (EPSG:3310) for metric accuracy, then project to template.
ca_3310   <- st_transform(ca_sf, "EPSG:3310")
ca_v_3310 <- vect(ca_3310)

# Coastline (Natural Earth)
coast_sf  <- rnaturalearth::ne_download(scale = 10, type = "coastline",
                                        category = "physical", returnclass = "sf")
coast_3310 <- st_transform(coast_sf, "EPSG:3310")
coast_v    <- vect(coast_3310)

# Metric template (~2 km) for distance calc, then project to 5'
tmpl_m <- rast(ext(ca_v_3310), crs = "EPSG:3310", resolution = 2000)
d_coast_m <- distance(tmpl_m, coast_v)                     # meters
d_coast_5 <- project(d_coast_m, template, method = "bilinear")
names(d_coast_5) <- "dist_coast_m"

# # NHDPlus flowlines (all stream flowlines in AOI), then distance in meters (US only)
#Last time I tried this was not working, doing google earth engine alternative below
# flowlines <- nhdplusTools::get_nhdplus(AOI = st_as_sf(ca_sf), realization = "flowline")
# flow_3310 <- st_transform(flowlines, "EPSG:3310")
# flow_v    <- vect(flow_3310)
# d_streams_m <- distance(tmpl_m, flow_v)
# d_streams_5 <- project(d_streams_m, template, method = "bilinear")
# names(d_streams_5) <- "dist_stream_m"
# 
# # Plot to verify
# plot(d_streams_5)
# plot(st_geometry(ca_outline), add = TRUE, border = "red", lwd = 2)

#### Distance to streams, including Mexico ####

  # Download Natural Earth rivers/streams dataset
rivers_sf <- rnaturalearth::ne_download(
    scale = 10, 
    type = "rivers_lake_centerlines",
    category = "physical", 
    returnclass = "sf"
  )

# Crop to your CA + 100km buffer extent (in WGS84)
rivers_ca <- st_crop(rivers_sf, st_bbox(ca_sf))

# Transform to California Albers for distance calculation
rivers_3310 <- st_transform(rivers_ca, "EPSG:3310")
rivers_v <- vect(rivers_3310)

# Calculate distance in meters (in Albers projection)
d_streams_m <- distance(tmpl_m, rivers_v)

# Project back to EPSG:4326
d_streams_wgs <- project(d_streams_m, "EPSG:4326", method = "bilinear")

# Crop to California + 100km buffer extent
d_streams_cropped <- terra::crop(d_streams_wgs, ca_v)

# Mask to the buffer area
d_streams_masked <- terra::mask(d_streams_cropped, ca_v, touches = FALSE)

# Resample to match your 5 arcmin template
d_streams_5 <- terra::resample(d_streams_masked, template, method = "bilinear")

# Name it
names(d_streams_5) <- "dist_stream_m"

# Plot to verify
plot(d_streams_5)
plot(st_geometry(ca_outline), add = TRUE, border = "red", lwd = 2)

# ---- Assemble predictor stack ----
core <- c(
  bio5_ca,                # WorldClim BIO1..BIO19 (5')
  elev_5m, slope_5m, aspect_5m, tri_5m,
  lc, hfp_5m, pop_5m,
  soil,
  d_coast_5, d_streams_5
)

# Clean names (WorldClim already bio1..bio19)
names(core) <- make.names(names(core), unique = TRUE)

# Crop/mask exactly to CA template extent
core <- crop(core, template)
core <- mask(core, ca_v, touches = TRUE)

# ---- Write stack + manifest ----
tif_out <- file.path(out_dir, "CA_SDM_predictors_5arcmin_epsg4326.tif")
writeRaster(core, tif_out, overwrite = TRUE,
            gdal = c("TILED=YES","COMPRESS=ZSTD","BIGTIFF=IF_SAFER"))

# Manifest CSV (layer, description, source)
desc <- tibble::tibble(
  layer = names(core),
  description = c(
    paste0("WorldClim v2.1 BIO", 1:19, " (5 arcmin)"),
    "Elevation (SRTM 30\" aggregated to 5')",
    "Slope (deg, aggregated to 5')",
    "Aspect (deg, aggregated to 5')",
    "Terrain Ruggedness Index (aggregated to 5')",
    paste("Land-cover fraction:", c("built","trees","shrubs","grassland",
                                    "cropland","bare","water","wetland")),
    "Landscape Intactness (2025, 1 km aggregated to 5')",
    "Human Footprint (2018, Venter et al. 2022)",
    "Population density (GPW v4, 2020, 5 arcmin)",
    "Soil sand % (0–30 cm, SoilGrids)",
    "Soil silt % (0–30 cm, SoilGrids)",
    "Soil clay % (0–30 cm, SoilGrids)",
    "Soil pH H₂O (0–30 cm, SoilGrids)",
    "Soil organic carbon (0–30 cm, SoilGrids)",
    "Soil bulk density (0–30 cm, SoilGrids)",
    "Distance to coast (m, EPSG:3310 → EPSG:4326)",
    "Distance to streams (m, Natural Earth rivers, EPSG:3310 → EPSG:4326)")[seq_len(nlyr(core))],
  source = c(
    rep("WorldClim 2.1 via geodata::worldclim_global", 19),
    "SRTM via geodata::elevation_30s (geomorphometry derived)",
    "Derived from SRTM",
    "Derived from SRTM",
    "Derived from SRTM",
    rep("ESA WorldCover fractions via geodata::landcover", 8),
    "Landscape Intactness (2025 Update, CBI GDB)",
    "Venter et al. 2022 Human Footprint Dataset",
    "GPW v4 Population Density (2020) via geodata::population",
    rep("SoilGrids (via geodata::soil_world)", 6),
    "Natural Earth Coastline (10 m, distance to coast)",
    "Natural Earth Rivers & Lake Centerlines (10 m, distance to streams)"
   
  )[seq_len(nlyr(core))]
)

readr::write_csv(desc, file.path(out_dir, "CA_SDM_predictors_manifest.csv"))

# Done: path summary
message("Wrote: ", tif_out)
message("Layers: ", nlyr(core))

# =========================
# California SDM - Enhanced Predictor Stack
# Complete End-to-End Script
# Adds: CWD, AET, VPD, Solar Radiation, LAI, NDVI, Fire History, Derived Variables
# =========================

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# STEP 1: RUN THIS SETUP CODE FIRST (once per R session):
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

library(reticulate)
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)

# Initialize Python EE
py_run_string("
import ee
ee.Initialize(project='species-distribution-modeling')
")

# Set environment variables
Sys.setenv(EARTHENGINE_PYTHON = py_config()$python)
Sys.setenv(RETICULATE_PYTHON = py_config()$python)

# Create session file
dir.create("~/.config/earthengine", recursive = TRUE, showWarnings = FALSE)
session_file <- path.expand("~/.config/earthengine/rgee_sessioninfo.txt")
cat(paste0("Python Path: ", py_config()$python, "\n"), file = session_file)

# Test GEE
library(rgee)
test <- ee$Image("USGS/SRTMGL1_003")
print(test$getInfo()$bands[[1]]$id)  # Should print "elevation"

message("✓ Setup complete! Now run the main script below.")

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# IF YOU SEE "elevation" ABOVE, YOU'RE READY TO RUN THE SCRIPT!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ---- Install & Load Packages ----
required <- c(
  "terra", "sf", "dplyr", "tidyr", "lubridate",
  "rgee", "reticulate",  # Google Earth Engine
  "MODISTools",          # MODIS data
  "httr", "jsonlite"     # For API calls
)

to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

lapply(required, library, character.only = TRUE)

setwd('/Users/chriscosma/Desktop/')

# ---- Paths & Setup ----
data_dir <- file.path(getwd(), "predictors_cache_enhanced")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

out_dir <- file.path(getwd(), "predictors_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load existing predictor stack to use as template
existing_stack_path <- file.path(out_dir, "CA_SDM_predictors_5arcmin_epsg4326.tif")

if (!file.exists(existing_stack_path)) {
  stop("Please run the original script first to create the base predictor stack!")
}

# Load template and CA boundary
core <- rast(existing_stack_path)
template <- core[[1]]  # 5 arcmin template in EPSG:4326

# ==== GRID LOCK (make EE outputs match the template exactly) ====
res_x  <- res(template)[1]
res_y  <- res(template)[2]
ext_t  <- ext(template)
xmin_t <- ext_t[1]; xmax_t <- ext_t[2]
ymin_t <- ext_t[3]; ymax_t <- ext_t[4]
ncol_t <- ncol(template); nrow_t <- nrow(template)

# EE region polygon as ring (xmin,ymin -> ... -> back to xmin,ymin)
region_coords <- list(list(
  c(xmin_t, ymin_t), c(xmin_t, ymax_t),
  c(xmax_t, ymax_t), c(xmax_t, ymin_t),
  c(xmin_t, ymin_t)
))

# (A) Push parameters to Python for py_run_string() calls
py$res_x         <- res_x
py$res_y         <- res_y
py$xmin_t        <- xmin_t
py$ymax_t        <- ymax_t
py$ncol_t        <- as.integer(ncol_t)
py$nrow_t        <- as.integer(nrow_t)
py$region_coords <- region_coords

# (B) R-side region geom for rgee calls
region_poly <- ee$Geometry$Polygon(region_coords, NULL, FALSE)

# California + 100km boundary
ca_sf <- rnaturalearth::ne_states(country = "United States of America",
                                  returnclass = "sf") |>
  filter(name == "California") |>
  st_make_valid() |>
  st_transform("EPSG:4326")

# Create 100km buffer in distance-preserving projection
ca_sf <- st_transform(ca_sf, "EPSG:3310")  # California Albers Equal Area
ca_sf <- st_buffer(ca_sf, dist = 100000)  # 100km = 100,000 meters

# Transform buffer back to EPSG:4326
ca_sf <- st_transform(ca_sf, "EPSG:4326")

# Convert both to SpatVector for terra operations
ca_v <- vect(ca_sf)
ca_bbox <- ext(ca_v)

# Define dynamic Earth Engine bounding box from full buffered CA extent
bbox_vals <- as.vector(st_bbox(ca_sf))
bbox <- ee$Geometry$Rectangle(
  list(bbox_vals[1], bbox_vals[2], bbox_vals[3], bbox_vals[4]),
  proj = "EPSG:4326",
  geodesic = FALSE
)
message("AOI extent (deg): ", paste(round(bbox_vals, 3), collapse = ", "))


message("Template loaded: ", res(template)[1], " degree resolution")
message("California extent: ", paste(as.vector(ca_bbox), collapse = ", "))

# =============================
# 1. TERRACLIMATE DATA (CWD, AET, VPD, Solar Radiation)
# Via Python-based download (bypasses rgee credential issues)
# =============================

message("\n=== Processing TerraClimate Data ===")

# Python-based download function
download_terraclim_py <- function(varname, band, year_start = 1991, year_end = 2020) {
  message("Downloading TerraClimate: ", varname, " (", year_start, "-", year_end, ")...")
  temp_file <- tempfile(fileext = ".tif")
  
  # provide variables to Python
  py$tc_band <- band
  py$ys      <- as.integer(year_start)
  py$ye      <- as.integer(year_end)
  
  py_run_string("
import ee

# Build region from R-provided coords
region = ee.Geometry.Polygon(region_coords, proj='EPSG:4326', geodesic=False)

# TerraClimate mean over years, select band
tc = (
    ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
    .filterDate(f'{ys}-01-01', f'{ye}-12-31')
    .select(tc_band)
    .mean()
)

# Reproject exactly onto the template grid (origin + cell size)
img = tc.reproject(
    crs='EPSG:4326',
    crsTransform=[res_x, 0, xmin_t, 0, -res_y, ymax_t]
).clip(region)

# Deterministic download (no 'scale'); use crsTransform + dimensions + region
url = img.getDownloadURL({
    'crs': 'EPSG:4326',
    'crsTransform': [res_x, 0, xmin_t, 0, -res_y, ymax_t],
    'dimensions': [ncol_t, nrow_t],
    'region': region,
    'format': 'GEO_TIFF'
})
")

# download and read
url <- py$url
download.file(url, temp_file, mode = "wb", quiet = FALSE)
r <- rast(temp_file)

# already on-grid; just crop/mask
r <- crop(r, ca_v) |> mask(ca_v, touches = TRUE)
names(r) <- paste0("terraclim_", tolower(gsub(" ", "_", varname)), "_", year_start, "_", year_end)

unlink(temp_file)
r
}

cwd_5m  <- download_terraclim_py("CWD",             "def",  1991, 2020)
aet_5m  <- download_terraclim_py("AET",             "aet",  1991, 2020)
vpd_5m  <- download_terraclim_py("VPD",             "vpd",  1991, 2020)
srad_5m <- download_terraclim_py("Solar Radiation", "srad", 1991, 2020)

terraclim_stack <- c(cwd_5m, aet_5m, vpd_5m, srad_5m)


###########################################################
### === Processing MODIS Data (Earth Engine, C7 updated) === ###
###########################################################

library(reticulate)
ee <- import("ee")
ee$Initialize()

tryCatch({
  message("=== Downloading MODIS NDVI and LAI (2018–2020, Collection 7) ===")
  
  # --- NDVI (MOD13Q1 C7) ---
  ndvi_ic <- ee$ImageCollection("MODIS/061/MOD13Q1")$
    filterDate("2018-01-01", "2020-12-31")$
    select("NDVI")
  
  ndvi_mean <- ndvi_ic$mean()$clip(bbox)
  ndvi_url <- ndvi_mean$getDownloadURL(list(
    scale = 4638,
    crs = "EPSG:4326",
    region = bbox,
    format = "GEO_TIFF"
  ))
  ndvi_path <- file.path(tempdir(), "modis_ndvi_2018_2020.tif")
  download.file(ndvi_url, ndvi_path, mode = "wb", quiet = FALSE)
  ndvi_rast <- rast(ndvi_path) |> resample(template, method = "bilinear") |> crop(ca_v) |> mask(ca_v, touches = TRUE)
  names(ndvi_rast) <- "modis_ndvi_2018_2020"
  
  # --- LAI (MOD15A2H C7) ---
  lai_ic <- ee$ImageCollection("MODIS/061/MOD15A2H")$
    filterDate("2018-01-01", "2020-12-31")$
    select("Lai_500m")
  
  lai_mean <- lai_ic$mean()$clip(bbox)
  lai_url <- lai_mean$getDownloadURL(list(
    scale = 4638,
    crs = "EPSG:4326",
    region = bbox,
    format = "GEO_TIFF"
  ))
  lai_path <- file.path(tempdir(), "modis_lai_2018_2020.tif")
  download.file(lai_url, lai_path, mode = "wb", quiet = FALSE)
  lai_rast  <- rast(lai_path)  |> resample(template, method = "bilinear") |> crop(ca_v) |> mask(ca_v, touches = TRUE)
  names(lai_rast) <- "modis_lai_2018_2020"
  
  modis_stack <- c(lai_rast, ndvi_rast)
  message("✓ MODIS data processed successfully (2 layers, Collection 7)")
  
}, error = function(e) {
  message("✗ Error with MODIS (Earth Engine): ", e$message)
  modis_stack <- NULL
})


###########################################################
### === Processing Fire History Data (MTBS) === ###
###########################################################

message("\n=== Processing Fire History Data (MTBS improved) ===")

tryCatch({
  mtbs_url <- "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/burned_area_extent_shapefile/mtbs_perimeter_data.zip"
  mtbs_zip <- file.path(data_dir, "mtbs_perimeters.zip")
  mtbs_dir <- file.path(data_dir, "mtbs_perimeters")
  
  if (!dir.exists(mtbs_dir)) {
    message("Downloading MTBS fire perimeters (~350 MB)...")
    download.file(mtbs_url, mtbs_zip, mode = "wb")
    unzip(mtbs_zip, exdir = mtbs_dir)
  }
  
  shp_files <- list.files(mtbs_dir, pattern = "\\.shp$", full.names = TRUE)
  if (length(shp_files) == 0) stop("No MTBS shapefiles found.")
  
  fires <- do.call(rbind, lapply(shp_files, terra::vect))
  fires <- project(fires, crs(template))
  fires_ca <- crop(fires, ca_v)
  
  # ---- Robust Year Field Detection ----
  yr_field <- intersect(
    c("MTBSFireYr", "Year", "YEAR", "Ig_Date"),
    names(fires_ca)
  )[1]
  if (is.na(yr_field)) stop("No recognizable fire year field found in MTBS data.")
  
  message("Using field for fire year: ", yr_field)
  
  # Extract year properly from SpatVector
  if (yr_field == "Ig_Date") {
    fires_ca$Year <- as.numeric(substr(fires_ca[[yr_field]][,1], 1, 4))
  } else {
    fires_ca$Year <- as.numeric(fires_ca[[yr_field]][,1])
  }
  
  # Filter using values() to access the data frame
  fire_df <- values(fires_ca)
  valid_idx <- which(!is.na(fire_df$Year) & fire_df$Year >= 1984 & fire_df$Year <= 2023)
  fires_ca <- fires_ca[valid_idx, ]
  
  if (nrow(fires_ca) == 0) stop("No fire records found within 1984-2023 range.")
  
  message("Found ", nrow(fires_ca), " fire perimeters in California (1984-2023)")
  
  # ---- Rasterize in Albers for Accuracy ----
  template_m <- project(template, "EPSG:3310")
  fires_alb  <- project(fires_ca, "EPSG:3310")
  
  fire_year  <- rasterize(fires_alb, template_m, field = "Year", fun = "max", background = NA)
  fire_freq  <- rasterize(fires_alb, template_m, fun = "count", background = 0)
  
  # ---- Project Back to EPSG:4326 ----
  fire_year_5 <- project(fire_year, template, method = "near")
  fire_freq_5 <- project(fire_freq, template, method = "near")
  
  # ---- Derived Variables ----
  time_since_fire <- 2024 - fire_year_5
  time_since_fire <- subst(time_since_fire, NA, 0)  # replace unburned with 0
  
  names(fire_freq_5)     <- "fire_frequency_1984_2023"
  names(time_since_fire) <- "time_since_fire_yrs"
  
  fire_stack <- c(fire_freq_5, time_since_fire) |>
    crop(ca_v) |>
    mask(ca_v, touches = TRUE)
  message("✓ Fire history data processed successfully (2 layers)")
  
}, error = function(e) {
  message("✗ Error processing fire data: ", e$message)
  fire_stack <- NULL
})

# =============================
# 4. CALCULATED VARIABLES
# Growing degree days, frost-free days, water balance approximations
# =============================

message("\n=== Calculating Derived Climate Variables ===")

# Identify WorldClim layers
bio1_idx  <- which(names(core) == "wc2.1_5m_bio_1")   # Mean annual temp
bio6_idx  <- which(names(core) == "wc2.1_5m_bio_6")   # Min temp coldest month
bio12_idx <- which(names(core) == "wc2.1_5m_bio_12")  # Annual precip
derived_stack <- NULL

if (length(bio1_idx) > 0 && length(bio6_idx) > 0 && length(bio12_idx) > 0) {
  bio1  <- core[[bio1_idx]]
  bio6  <- core[[bio6_idx]]
  bio12 <- core[[bio12_idx]]
  
  # --- Growing degree days (base 10°C & 5°C) ---
  gdd_base10 <- app(bio1, fun = function(x) pmax(x - 10, 0) * 365)
  names(gdd_base10) <- "gdd_base10_annual"
  
  gdd_base5 <- app(bio1, fun = function(x) pmax(x - 5, 0) * 365)
  names(gdd_base5) <- "gdd_base5_annual"
  
  # --- Frost-free days estimate ---
  frost_free <- app(bio6, fun = function(x) 365 * pmax(x + 10, 0) / 30)
  frost_free <- clamp(frost_free, lower = 0, upper = 365)
  names(frost_free) <- "frost_free_days_est"
  
  # --- Simplified water balance (fallback only if TerraClimate failed) ---
  if (is.null(terraclim_stack)) {
    message("Calculating simplified water balance from WorldClim...")
    
    # Thornthwaite PET approximation (mm)
    pet_mm <- app(bio1, fun = function(t) {
      t_adj <- pmax(t, 0)
      16 * 12 * (t_adj / 5)^1.5
    })
    names(pet_mm) <- "pet_thornthwaite_mm"
    
    # Simplified AET (limited by precip + storage)
    aet_mm <- overlay(bio12, pet_mm, fun = function(p, pet) pmin(p + 100, pet))
    names(aet_mm) <- "aet_simplified_mm"
    
    # Climatic Water Deficit
    cwd_mm <- overlay(pet_mm, aet_mm, fun = function(pet, aet) pmax(pet - aet, 0))
    names(cwd_mm) <- "cwd_simplified_mm"
    
    derived_stack <- c(gdd_base10, gdd_base5, frost_free, pet_mm, aet_mm, cwd_mm)
  } else {
    derived_stack <- c(gdd_base10, gdd_base5, frost_free)
  }
  
  message("✓ Derived variables calculated (", nlyr(derived_stack), " layers)")
  
} else {
  message("✗ Cannot calculate derived variables – missing temperature data")
}


###########################################################
### === Processing Fire-CCI v5.1 Data (ESA via GEE) === ###
###########################################################

message("\n=== Processing Fire-CCI v5.1 Data (ESA) ===")

library(reticulate)
ee <- import("ee")
ee$Initialize()

tryCatch({
  message("=== Downloading Fire-CCI v5.1 (2001–2020) ===")
  
  # --- Fire-CCI v5.1 burned area product ---
  # Collection: ESA/CCI/FireCCI/5_1
  # Band: BurnDate (day of year when burned, or 0 if not burned)
  fire_ic <- ee$ImageCollection("ESA/CCI/FireCCI/5_1")$
    filterDate("2001-01-01", "2020-12-31")$
    filterBounds(bbox)
  
  # Select only the BurnDate band
  fire_ic <- fire_ic$select("BurnDate")
  
  # --- Calculate fire frequency (number of times pixel burned) ---
  # Create binary mask (1 = burned, 0 = not burned)
  binary_fires <- fire_ic$map(function(img) {
    img$gt(0)  # BurnDate > 0 means burned
  })
  
  fire_frequency <- binary_fires$sum()$clip(bbox)
  
  # Align the summed binary fires onto the template grid before download
  fire_frequency_on_grid <- fire_frequency$reproject(
    crs = "EPSG:4326",
    crsTransform = list(res_x, 0, xmin_t, 0, -res_y, ymax_t)
  )$clip(region_poly)
  
  fire_freq_url <- fire_frequency_on_grid$getDownloadURL(list(
    crs = "EPSG:4326",
    crsTransform = list(res_x, 0, xmin_t, 0, -res_y, ymax_t),
    dimensions = list(as.integer(ncol_t), as.integer(nrow_t)),
    region = region_poly,
    format = "GEO_TIFF"
  ))
  
  fire_freq_path <- file.path(tempdir(), "firecci_frequency_2001_2020.tif")
  download.file(fire_freq_url, fire_freq_path, mode = "wb", quiet = FALSE)
  fire_freq_rast <- rast(fire_freq_path) |>
    crop(ca_v) |>
    mask(ca_v, touches = TRUE)
  names(fire_freq_rast) <- "firecci_frequency_2001_2020"
  
  
  # --- Most recent burn year (time since fire) ---
  # Get the maximum (most recent) BurnDate, then extract year
  most_recent_fire <- fire_ic$max()$clip(bbox)
  
  # Convert day-of-year to year (requires the year property from metadata)
  # Alternative: use a reducer to get the latest image's year
  latest_fire_year <- fire_ic$
    map(function(img) {
      year <- ee$Number(ee$Date(img$get("system:time_start"))$get("year"))
      # Create constant year image, cast to Int16 so all have same type
      year_img <- ee$Image$constant(year)$toInt16()$rename("year")
      return(year_img$updateMask(img$gt(0)))  # mask to burned pixels
    })$
    max()$
    clip(bbox)
  
  
  latest_fire_on_grid <- latest_fire_year$reproject(
    crs = "EPSG:4326",
    crsTransform = list(res_x, 0, xmin_t, 0, -res_y, ymax_t)
  )$clip(region_poly)
  
  latest_fire_url <- latest_fire_on_grid$getDownloadURL(list(
    crs = "EPSG:4326",
    crsTransform = list(res_x, 0, xmin_t, 0, -res_y, ymax_t),
    dimensions = list(as.integer(ncol_t), as.integer(nrow_t)),
    region = region_poly,
    format = "GEO_TIFF"
  ))
  
  latest_fire_path <- file.path(tempdir(), "firecci_latest_year_2001_2020.tif")
  download.file(latest_fire_url, latest_fire_path, mode = "wb", quiet = FALSE)
  latest_fire_rast <- rast(latest_fire_path) |>
    crop(ca_v) |>
    mask(ca_v, touches = TRUE)
  
  # Calculate time since fire (2024 - last fire year)
  # --- Fix: compute time since fire correctly ---
  # Replace 0 (unburned) with NA BEFORE subtraction
  latest_fire_rast[latest_fire_rast == 0] <- NA
  
  # Compute years since last fire
  time_since_fire <- 2024 - latest_fire_rast
  names(time_since_fire) <- "firecci_time_since_fire_yrs"
  
  # Optionally, for modeling consistency, fill unburned cells with 0 (not 100)
  time_since_fire <- subst(time_since_fire, NA, 0)
  
  # --- Total burned area fraction (optional) ---
  # Calculate proportion of years pixel was burned
  fire_proportion <- fire_freq_rast / 20
  names(fire_proportion) <- "firecci_burn_proportion_2001_2020"
  
  # Combine fire layers
  firecci_stack <- c(fire_freq_rast, time_since_fire, fire_proportion)
  
  message("✓ Fire-CCI v5.1 data processed successfully (3 layers)")
  message("  - Fire frequency: ", minmax(fire_freq_rast)[1], " to ", minmax(fire_freq_rast)[2], " fires")
  message("  - Time since fire: ", minmax(time_since_fire)[1], " to ", minmax(time_since_fire)[2], " years")
  
}, error = function(e) {
  message("✗ Error with Fire-CCI: ", e$message)
  message("Skipping Fire-CCI data.")
  firecci_stack <- NULL
})

###########################################################
### === Distance to Water (JRC Global Surface Water) === ###
###########################################################

message("\n=== Processing Distance to Water (JRC Global Surface Water) ===")

library(reticulate)
ee <- import("ee")
ee$Initialize()

tryCatch({
  message("=== Downloading water bodies via Google Earth Engine ===")
  
  # --- JRC Global Surface Water (includes rivers, lakes, ponds) ---
  # Use "occurrence" band: % of time water was present (1984-2021)
  # Threshold at >50% occurrence to define persistent water
  jrc_water <- ee$Image("JRC/GSW1_4/GlobalSurfaceWater")$
    select("occurrence")$
    gt(50)$  # Water present >50% of the time
    clip(bbox)
  
  # Download water mask
  water_url <- jrc_water$getDownloadURL(list(
    scale = 1000,  # 1km for faster processing, will resample later
    crs = "EPSG:4326",
    region = bbox,
    format = "GEO_TIFF"
  ))
  
  water_path <- file.path(tempdir(), "jrc_water_mask.tif")
  download.file(water_url, water_path, mode = "wb", quiet = FALSE)
  water_mask <- rast(water_path)
  
  message("Calculating distance to water in California Albers projection...")
  
  # Convert to California Albers for accurate distance calculation
  water_mask_alb <- project(water_mask, "EPSG:3310", method = "near")
  
  # Create points from water pixels (where value = 1)
  water_mask_alb[water_mask_alb == 0] <- NA
  water_pts <- as.points(water_mask_alb, values = TRUE, na.rm = TRUE)
  
  # Calculate distance to nearest water point
  template_alb <- project(template, "EPSG:3310")
  d_water_m <- distance(template_alb, water_pts)
  
  # Project back to EPSG:4326 and align to template
  d_water_5 <- project(d_water_m, template, method = "bilinear")
  d_water_5 <- crop(d_water_5, ca_v)
  d_water_5 <- mask(d_water_5, ca_v, touches = FALSE)
  
  names(d_water_5) <- "dist_water_m_jrc"
  
  message("✓ Distance to water calculated successfully")
  message("  - Min distance: ", round(minmax(d_water_5)[1]), " m")
  message("  - Max distance: ", round(minmax(d_water_5)[2]), " m")
  
}, error = function(e) {
  message("✗ Error with JRC water distance: ", e$message)
  d_water_5 <- NULL
})

###########################################################
### === Distance to Streams Only (HydroSHEDS) === ###
###########################################################

message("\n=== Processing Distance to Streams (HydroSHEDS) ===")

tryCatch({
  message("=== Downloading stream network via Google Earth Engine ===")
  
  # HydroSHEDS river network (15s resolution ~500m)
  # Use only rivers with upstream area > 100 km² for major streams
  rivers <- ee$Image("WWF/HydroSHEDS/15ACC")$  # Flow accumulation
    gt(200)$  # Threshold for stream definition (adjust as needed)
    clip(bbox)
  
  # Download stream mask
  rivers_url <- rivers$getDownloadURL(list(
    scale = 1000,
    crs = "EPSG:4326",
    region = bbox,
    format = "GEO_TIFF"
  ))
  
  rivers_path <- file.path(tempdir(), "hydrosheds_rivers.tif")
  download.file(rivers_url, rivers_path, mode = "wb", quiet = FALSE)
  rivers_mask <- rast(rivers_path)
  
  message("Calculating distance to streams in California Albers projection...")
  
  # Same distance calculation as above
  rivers_mask_alb <- project(rivers_mask, "EPSG:3310", method = "near")
  rivers_mask_alb[rivers_mask_alb == 0] <- NA
  rivers_pts <- as.points(rivers_mask_alb, values = TRUE, na.rm = TRUE)
  
  template_alb <- project(template, "EPSG:3310")
  d_streams_m <- distance(template_alb, rivers_pts)
  
  d_streams_5 <- project(d_streams_m, template, method = "bilinear")
  d_streams_5 <- crop(d_streams_5, ca_v)
  d_streams_5 <- mask(d_streams_5, ca_v, touches = FALSE)
  
  names(d_streams_5) <- "dist_stream_m_hydrosheds"
  
  message("✓ Distance to streams calculated successfully")
  
}, error = function(e) {
  message("✗ Error with HydroSHEDS stream distance: ", e$message)
  d_streams_5 <- NULL
})

# =============================
# 5. COMBINE ALL NEW LAYERS
# =============================

message("\n=== Combining Enhanced Predictor Stack ===")

# Collect all new layers
new_layers <- list(
  terraclim_stack,
  modis_stack,
  fire_stack,
  firecci_stack,   # ADD FireCCI stack
  d_water_5,       # ADD JRC distance to water
  d_streams_5,     # ADD HydroSHEDS distance to streams
  derived_stack
)


# Remove NULL elements
new_layers <- new_layers[!sapply(new_layers, is.null)]

if (length(new_layers) > 0) {
  # Combine new layers
  enhanced_predictors <- do.call(c, new_layers)
  ca_v <- project(ca_v, crs(template))
  enhanced_predictors <- crop(enhanced_predictors, ca_v)
  enhanced_predictors <- mask(enhanced_predictors, ca_v, touches = TRUE)
  
  # Combine with existing stack
  full_stack <- c(core, enhanced_predictors)
  names(full_stack) <- make.names(names(full_stack), unique = TRUE)
  
  message("\n=== Summary ===")
  message("Original layers: ", nlyr(core))
  message("New layers added: ", nlyr(enhanced_predictors))
  message("Total layers: ", nlyr(full_stack))
  
  # ---- Write Enhanced Stack ----
  tif_out <- file.path(out_dir, "CA_SDM_predictors_5arcmin_ENHANCED.tif")
  message("\nWriting enhanced stack to disk...")
  writeRaster(full_stack, tif_out, overwrite = TRUE,
              gdal = c("TILED=YES", "COMPRESS=ZSTD", "BIGTIFF=IF_SAFER"))
  
  # ---- Create Manifest ----
  manifest <- data.frame(
    layer = names(enhanced_predictors),
    description = c(
      if (!is.null(terraclim_stack)) c(
        "Climatic water deficit - TerraClimate 1991–2020 mean (mm)",
        "Actual evapotranspiration - TerraClimate 1991–2020 mean (mm)",
        "Vapor pressure deficit - TerraClimate 1991–2020 mean (kPa)",
        "Solar radiation - TerraClimate 1991–2020 mean (W/m²)"
      ),
      if (!is.null(modis_stack)) c(
        "Leaf Area Index (LAI) - MODIS 2018–2020 mean",
        "Normalized Difference Vegetation Index (NDVI) - MODIS 2018–2020 mean"
      ),
      if (!is.null(fire_stack)) c(
        "Fire frequency (1984–2023) - MTBS perimeters",
        "Time since last fire (years since 2024) - MTBS perimeters"
      ),
      if (!is.null(firecci_stack)) c(
        "Fire frequency (2001–2020) - ESA FireCCI v5.1",
        "Time since last fire (years since 2024) - ESA FireCCI v5.1",
        "Burned fraction (proportion of years burned) - ESA FireCCI v5.1"
      ),
      if (!is.null(d_water_5)) "Distance to persistent surface water (m, JRC GSW 1984–2021)",
      if (!is.null(d_streams_5)) "Distance to major streams (m, HydroSHEDS ≥100 km²)",
      if (!is.null(derived_stack)) {
        desc <- c(
          "Growing degree days base 10°C (annual)",
          "Growing degree days base 5°C (annual)",
          "Frost-free days (estimated)"
        )
        if (is.null(terraclim_stack)) {
          desc <- c(desc,
                    "Potential evapotranspiration - Thornthwaite (mm)",
                    "Actual evapotranspiration - simplified (mm)",
                    "Climatic water deficit - simplified (mm)")
        }
        desc
      }
    )[1:nlyr(enhanced_predictors)],
    source = c(
      if (!is.null(terraclim_stack)) rep("TerraClimate via GEE", 4),
      if (!is.null(modis_stack)) rep("MODIS (NASA) via GEE", 2),
      if (!is.null(fire_stack)) rep("MTBS (USGS, 1984–2023)", 2),
      if (!is.null(firecci_stack)) rep("ESA FireCCI v5.1 (2001–2020)", 3),
      if (!is.null(d_water_5)) "JRC Global Surface Water (1984–2021)",
      if (!is.null(d_streams_5)) "HydroSHEDS 15s (WWF)",
      if (!is.null(derived_stack)) {
        src <- rep("Calculated from WorldClim", 3)
        if (is.null(terraclim_stack)) {
          src <- c(src, rep("Calculated from WorldClim", 3))
        }
        src
      }
    )[1:nlyr(enhanced_predictors)]
  )
  
  write.csv(manifest, 
            file.path(out_dir, "CA_SDM_predictors_manifest_ENHANCED.csv"),
            row.names = FALSE)
  
  # ---- Preview Plots ----
  message("Creating preview plots...")
  png(file.path(out_dir, "enhanced_predictors_preview.png"),
      width = 2400, height = 1800, res = 150)
  
  n_plots <- nlyr(enhanced_predictors)
  par(mfrow = c(ceiling(n_plots/4), 4), mar = c(2, 2, 3, 2))
  
  for (i in 1:n_plots) {
    plot(enhanced_predictors[[i]], 
         main = names(enhanced_predictors)[i],
         col = hcl.colors(64, "Spectral", rev = TRUE),
         axes = FALSE)
    plot(ca_v, add = TRUE, border = "gray30", lwd = 0.5)
  }
  
  dev.off()
  
  message("\n✓✓✓ COMPLETE! ✓✓✓")
  message("Enhanced stack saved: ", tif_out)
  message("Manifest saved: CA_SDM_predictors_manifest_ENHANCED.csv")
  message("Preview plots saved: enhanced_predictors_preview.png")
  message("\nNew predictors summary:")
  print(manifest[, c("layer", "source")])
  
} else {
  message("\n✗ WARNING: No new layers were successfully created.")
  message("Check error messages above for details.")
}

message("\n=== Script Complete ===")

message("\nNA diagnostics after masking:")
na_prop <- terra::global(full_stack, fun = function(x, na.rm=FALSE) sum(is.na(x))/length(x))
print(round(na_prop, 3))

#### Final fixes ####

# --- 1) Read the stack you already wrote ---
out_dir <- file.path(getwd(), "predictors_out")
enhanced_path <- file.path(out_dir, "CA_SDM_predictors_5arcmin_ENHANCED.tif")
base_path     <- file.path(out_dir, "CA_SDM_predictors_5arcmin_epsg4326.tif")
env_stack <- if (file.exists(enhanced_path)) terra::rast(enhanced_path) else terra::rast(base_path)

# --- 2) Mask everything to the non-NA footprint of WorldClim BIO1 ---
# robustly find the BIO1 layer name
bio1_name <- names(env_stack)[grep("bio[_\\.]*1$", names(env_stack))]
stopifnot(length(bio1_name) == 1)
bio1 <- env_stack[[bio1_name]]

# align extents just in case, then mask
env_stack_masked <- terra::mask(env_stack, bio1)  # keeps cells where BIO1 is not NA

# --- 3) Gap-fill ONLY within the BIO1 footprint ---
# Best-practice, small-gap local interpolation:
#   pass 1: 3x3 median (robust to outliers)
#   pass 2: 5x5 mean (smooths residual speckle)
fill_one <- function(x) {
  # operate only where x is NA but BIO1 is not NA
  # (terra::focal with na.policy="only" computes output only for NA cells)
  x1 <- terra::cover(x, terra::focal(x, w = 3, fun = median,
                                     na.rm = TRUE, na.policy = "only"))
  x2 <- terra::cover(x1, terra::focal(x1, w = 5, fun = mean,
                                      na.rm = TRUE, na.policy = "only"))
  x2
}

env_filled <- terra::rast(lapply(seq_len(terra::nlyr(env_stack_masked)),
                                 function(i) fill_one(env_stack_masked[[i]])))
names(env_filled) <- names(env_stack_masked)

# Ensure we did not fill outside BIO1’s valid area
env_stack_final <- terra::mask(env_filled, bio1)

# Layers to drop
drop_layers <- c(
  "firecci_frequency_2001_2020",
  "firecci_time_since_fire_yrs",
  "firecci_burn_proportion_2001_2020",
  "time_since_fire_yrs",
  "lc_water",
  "dist_stream_m",
  "dist_stream_m_hydrosheds"
)

# Remove any that exist in the stack
keep_layers <- setdiff(names(env_stack_final), drop_layers)
env_stack_final <- env_stack_final[[keep_layers]]

# --- 4) Write final stack ---
tif_final <- file.path(out_dir, "env_stack_final.tif")
terra::writeRaster(env_stack_final, tif_final, overwrite = TRUE,
                   gdal = c("TILED=YES","COMPRESS=ZSTD","BIGTIFF=IF_SAFER"))
message("✓ Wrote: ", tif_final)

names(env_stack_final)

###########################################################
### === Quick Look at Enhanced Environmental Stack === ###
###########################################################

# Path to enhanced or base stack

env_stack = rast(file.path(out_dir, "env_stack_final.tif"))

stopifnot(inherits(env_stack, "SpatRaster"))

# Output folder for quicklook images
out_dir_imgs <- file.path(out_dir, "quicklook_env_layers")
dir.create(out_dir_imgs, showWarnings = FALSE, recursive = TRUE)

# Simple palette
pal <- hcl.colors(64, "YlGnBu", rev = TRUE)

# Check for California boundary vector
has_ca <- exists("ca_v") && inherits(ca_v, "SpatVector")

n   <- nlyr(env_stack)
nms <- names(env_stack)

message("Generating ", n, " quicklook PNGs...")
pb <- txtProgressBar(min = 0, max = n, style = 3)

for (i in seq_len(n)) {
  nm  <- make.names(nms[i])
  lyr <- env_stack[[i]]
  f   <- file.path(out_dir_imgs, sprintf("%03d_%s.png", i, nm))
  
  # Sanitize filenames
  f <- gsub("[^A-Za-z0-9_\\-\\.]", "_", f)
  
  png(f, width = 1600, height = 1200, res = 150)
  par(mar = c(2.2, 2.2, 2.5, 3))
  
  # Fast plot without unsupported args
  terra::plot(
    lyr,
    col    = pal,
    colNA  = "red",  
    main   = nm,
    axes   = FALSE,
    legend = TRUE,
    plg    = list(bg = "white")  # keep legend readable; no frame option
  )
  
  if (has_ca) terra::plot(ca_v, add = TRUE, col = NA, border = "grey20", lwd = 1)
  
  dev.off()
  setTxtProgressBar(pb, i)
}
close(pb)


message("\n✓ Wrote ", n, " quicklook PNGs to: ", out_dir_imgs)

names(env_stack_final)



