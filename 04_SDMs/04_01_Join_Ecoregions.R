#### Ecoregions (create once) ####

# --- Jepson ecoregions (California only) ---
jepson_ecoregions <- st_read(
  "Data_Clean/SDMs/SDM_Inputs/jepson_ecoregions/Jepson_regions_archive/CA-SUBDIV-v2-5-3_07012022-4326.shp",
  quiet = TRUE
)
jepson_proj <- st_transform(jepson_ecoregions, ca_albers)

# --- California boundary ---
ca_boundary <- geodata::gadm(country = "USA", level = 1, path = data_dir)
ca_boundary <- ca_boundary[ca_boundary$NAME_1 == "California", ]
ca_boundary_proj <- st_transform(st_as_sf(ca_boundary), ca_albers)

# --- EPA Level III ecoregions (North America) ---
epa_ecoregions <- st_read(
  "Data_Clean/SDMs/SDM_Inputs/environmental/NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp",
  quiet = TRUE
)
epa_proj <- st_transform(epa_ecoregions, ca_albers)

# --- 100 km buffer outside California ---
ca_buffer_100km <- st_buffer(ca_boundary_proj, dist = 100000)
buffer_zone_only <- st_difference(ca_buffer_100km, ca_boundary_proj)

# --- Clip EPA ecoregions to just the buffer zone ---
epa_buffer <- st_intersection(epa_proj, buffer_zone_only)

# --- Standardize attribute names before merging ---
jepson_proj_slim <- jepson_proj %>%
  dplyr::select(JEPCODE, geometry) %>%
  rename(ecoregion_code = JEPCODE) %>%
  mutate(source = "Jepson")

epa_buffer_slim <- epa_buffer %>%
  dplyr::select(NA_L3CODE, geometry) %>%
  rename(ecoregion_code = NA_L3CODE) %>%
  mutate(source = "EPA")

# --- Combine into one object ---
jepson_proj_combined <- dplyr::bind_rows(jepson_proj_slim, epa_buffer_slim)

# --- Validate geometries and fix issues if needed ---
jepson_proj_combined <- st_make_valid(jepson_proj_combined)

# --- Plot to verify ---
# Separate the two sources for coloring
jepson_only <- jepson_proj_combined %>% filter(source == "Jepson")
epa_only    <- jepson_proj_combined %>% filter(source == "EPA")

# Plot
# Compute combined extent (bounding box)
combined_bbox <- st_bbox(rbind(jepson_only, epa_only))

# Plot empty frame first to ensure full coverage
plot(0, 0, type = "n",
     xlim = c(combined_bbox["xmin"], combined_bbox["xmax"]),
     ylim = c(combined_bbox["ymin"], combined_bbox["ymax"]),
     asp = 1,
     xlab = "", ylab = "",
     main = "California + Buffer Ecoregions")

# Then add the layers
plot(st_geometry(jepson_only), col = "lightgreen", border = "grey40", add = TRUE)
plot(st_geometry(epa_only), col = "lightblue", border = "grey40", add = TRUE)

#save
# 1. Make sure the directory exists
dir.create("Data_Clean/SDMs/SDM_Inputs/environmental/combined_ecoregions",
           showWarnings = FALSE, recursive = TRUE)

# 2. Use a colon-free directory name (recommended)
out_path <- "Data_Clean/SDMs/SDM_Inputs/environmental/combined_ecoregions/jepson_epa_ecoregions_combined.gpkg"

# 3. Explicitly set the driver
st_write(jepson_proj_combined, out_path, driver = "GPKG", delete_dsn = TRUE)
message("✅ Saved combined ecoregion GeoPackage to: ", out_path)