#### Step one: Go to https://www.calflora.org/entry/observ.html, download occurrences for somewhere (Tools -> Download at the top), and extract the wkt value (API key) from the download link before it fully downloads (works in chrome, not firefox)

mtk_token = ""

#### Load required libraries ####
library(sf)
library(dplyr)
library(httr)
library(stringr)
library(tigris)
library(purrr)
library(readr)
library(maps)
library(ggplot2)

#### Functions ####
# Create WKT polygon from bounding box
create_wkt_polygon <- function(bbox) {
  # bbox should be c(xmin, ymin, xmax, ymax)
  # CalFlora expects: lon lat,lon lat,... format (space-separated, not comma between lon/lat)
  wkt <- sprintf("%.5f+%.5f,%.5f+%.5f,%.5f+%.5f,%.5f+%.5f,%.5f+%.5f",
                 bbox[1], bbox[4],  # top-left (xmin, ymax)
                 bbox[3], bbox[4],  # top-right (xmax, ymax)
                 bbox[3], bbox[2],  # bottom-right (xmax, ymin)
                 bbox[1], bbox[2],  # bottom-left (xmin, ymin)
                 bbox[1], bbox[4])  # close polygon (xmin, ymax)
  return(wkt)
}

# Download CalFlora data
download_calflora_data <- function(wkt, output_file, base_url, max_retries = 3) {
  # Construct the full URL with the WKT parameter
  url <- paste0(base_url, "&wkt=", URLencode(wkt, reserved = TRUE))
  
  # Try downloading with retries
  for (attempt in 1:max_retries) {
    tryCatch({
      response <- GET(url, timeout(120))
      
      if (status_code(response) == 200) {
        writeBin(content(response, "raw"), output_file)
        cat(sprintf("✓ Downloaded: %s (Attempt %d)\n", output_file, attempt))
        return(TRUE)
      } else {
        cat(sprintf("✗ HTTP %d for %s (Attempt %d)\n", 
                    status_code(response), output_file, attempt))
      }
    }, error = function(e) {
      cat(sprintf("✗ Error for %s (Attempt %d): %s\n", output_file, attempt, e$message))
    })
    
    if (attempt < max_retries) Sys.sleep(2)
  }
  
  return(FALSE)
}

# Main function to divide state and download
download_calflora_grid <- function(n_chunks = 500, output_dir = "calflora_data") {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## Create accurate California + 100km buffer using CA Albers
  cat("Creating California boundary with 100km buffer...\n")
  
  # Get California state boundary
  california_sf <- states(cb = TRUE, resolution = "500k") %>%
    filter(STUSPS == "CA") %>%
    st_make_valid()
  
  # Project to CA Albers Equal Area (EPSG:3310) for accurate distance
  california_proj <- st_transform(california_sf, crs = 3310)
  
  # Create 100km buffer (100,000 meters)
  california_buffer <- st_buffer(california_proj, dist = 100000)
  
  # Transform back to WGS84 for CalFlora API
  california_buffer_wgs84 <- st_transform(california_buffer, crs = 4326)
  
  # Get bounding box of the buffered area
  ca_bbox <- st_bbox(california_buffer_wgs84)
  
  cat(sprintf("Buffered area bbox: xmin=%.4f, ymin=%.4f, xmax=%.4f, ymax=%.4f\n",
              ca_bbox["xmin"], ca_bbox["ymin"], ca_bbox["xmax"], ca_bbox["ymax"]))
  
  ## Calculate grid dimensions
  total_area <- (ca_bbox["xmax"] - ca_bbox["xmin"]) * (ca_bbox["ymax"] - ca_bbox["ymin"])
  cell_area <- total_area / n_chunks
  cell_side <- sqrt(cell_area)
  
  n_cols <- ceiling((ca_bbox["xmax"] - ca_bbox["xmin"]) / cell_side)
  n_rows <- ceiling((ca_bbox["ymax"] - ca_bbox["ymin"]) / cell_side)
  
  cat(sprintf("Creating %d x %d grid (%d total cells)\n", n_cols, n_rows, n_cols * n_rows))
  
  # Base URL (everything except  WKT parameter)
  base_url <- paste0("https://calflora.org/app/download?cols=ID,Plant,Photo,Observer,Date,Source,County,Location+Description,Access,Accuracy:+Square+Meters,Citation,Calrecnum,Common+Name,Cover,Dataset,Date/Time,Datum,Detail+URL,Distribution,Distribution+Code,Elevation+(m),Entered+by,Error+Radius+(m),From,Genus,Geometry+Id,Gross+Area+Count,Gross+Area+Units,GrossArea,Group+Number,Habitat,Herbarium,Index+Date,Infested+Area+Count,Infested+Area+Units,InfestedArea,Large+Photo,Latitude,Location+Quality,Longitude,Management+Status,Modify+Date,Most+Recent,National+Ownership,National+Ownership+Code,Native+Status,Natural+Status,Notes,Number+of+Plants,Observation+Date,Organization,Percent+Cover,Percent+Cover+Midpoint,Phenology,Phenology+Code,Photo+URL,Region,Root,Species,Subspecies,Survey+ID,Taxon&format=CSV&mtk=", mtk_token, "&cell=t&georef=xc&georeferenced=t&ostatus=a&pcount=1&doc_type=rs&natural_status=w&nstatus=2&addnloc=L&cch=t&cnabh=t&wint=r")
  
  ## Generate grid cells and download
  successful_downloads <- 0
  failed_downloads <- 0
  
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      # Calculate cell bounds
      xmin <- ca_bbox["xmin"] + (j - 1) * cell_side
      xmax <- min(ca_bbox["xmin"] + j * cell_side, ca_bbox["xmax"])
      ymin <- ca_bbox["ymin"] + (i - 1) * cell_side
      ymax <- min(ca_bbox["ymin"] + i * cell_side, ca_bbox["ymax"])
      
      cell_bbox <- c(xmin, ymin, xmax, ymax)
      wkt <- create_wkt_polygon(cell_bbox)
      
      # Create filename
      output_file <- file.path(output_dir, sprintf("calflora_row%03d_col%03d.csv", i, j))
      
      # Download data
      cat(sprintf("[%d/%d] Processing cell (%d, %d)...\n", 
                  (i-1)*n_cols + j, n_cols*n_rows, i, j))
      
      if (download_calflora_data(wkt, output_file, base_url)) {
        successful_downloads <- successful_downloads + 1
      } else {
        failed_downloads <- failed_downloads + 1
      }
      
      # Be respectful to the server
      Sys.sleep(1)
    }
  }
  
  cat(sprintf("\n=== Download Complete ===\n"))
  cat(sprintf("Successful: %d\n", successful_downloads))
  cat(sprintf("Failed: %d\n", failed_downloads))
  cat(sprintf("Output directory: %s\n", output_dir))
}

# Check grid coverage
check_calflora_grid <- function(data_dir = "Data_Raw/Occurrences/calflora_data",
                                problem_cells_file_dir = "Temp/",
                                n_chunks = 500) {
  
  # Get all downloaded files
  files <- list.files(data_dir, pattern = "calflora_row.*\\.csv$", full.names = TRUE)
  
  cat(sprintf("Found %d files in directory\n", length(files)))
  
  # Extract row and column numbers from filenames
  file_info <- data.frame(
    file = files,
    basename = basename(files)
  ) %>%
    mutate(
      row = as.integer(str_extract(basename, "(?<=row)\\d+")),
      col = as.integer(str_extract(basename, "(?<=col)\\d+"))
    )
  
  # Check file sizes and record counts
  file_info <- file_info %>%
    mutate(
      size_bytes = file.size(file),
      exists = file.exists(file)
    )
  
  # Read each file to count records
  cat("Checking file contents...\n")
  file_info <- file_info %>%
    mutate(
      n_records = sapply(file, function(f) {
        tryCatch({
          df <- read_csv(f, col_types = cols(.default = "c"), show_col_types = FALSE)
          nrow(df)
        }, error = function(e) {
          NA_integer_
        })
      })
    )
  
  # Reconstruct the grid dimensions
  california_sf <- states(cb = TRUE, resolution = "500k") %>%
    filter(STUSPS == "CA") %>%
    st_make_valid()
  
  california_proj <- st_transform(california_sf, crs = 3310)
  california_buffer <- st_buffer(california_proj, dist = 100000)
  california_buffer_wgs84 <- st_transform(california_buffer, crs = 4326)
  ca_bbox <- st_bbox(california_buffer_wgs84)
  
  total_area <- (ca_bbox["xmax"] - ca_bbox["xmin"]) * (ca_bbox["ymax"] - ca_bbox["ymin"])
  cell_area <- total_area / n_chunks
  cell_side <- sqrt(cell_area)
  
  n_cols <- ceiling((ca_bbox["xmax"] - ca_bbox["xmin"]) / cell_side)
  n_rows <- ceiling((ca_bbox["ymax"] - ca_bbox["ymin"]) / cell_side)
  
  # Create complete grid
  complete_grid <- expand.grid(
    row = 1:n_rows,
    col = 1:n_cols
  )
  
  # Join with actual files
  grid_status <- complete_grid %>%
    left_join(file_info, by = c("row", "col")) %>%
    mutate(
      status = case_when(
        is.na(exists) ~ "missing",
        size_bytes < 100 ~ "empty_or_corrupt",
        is.na(n_records) ~ "read_error",
        n_records == 0 ~ "no_records",
        TRUE ~ "ok"
      )
    )
  
  # Summary statistics
  cat("\n=== Grid Coverage Summary ===\n")
  print(table(grid_status$status))
  
  # Identify problem cells
  problem_cells <- grid_status %>%
    filter(status != "ok") %>%
    arrange(row, col)
  
  if (nrow(problem_cells) > 0) {
    cat("\n=== Problem Cells ===\n")
    print(problem_cells %>% select(row, col, status, size_bytes, n_records))
    
    # Save to CSV for easy reprocessing
    write.csv(problem_cells, paste0(problem_cells_file_dir, "calflora_problem_cells.csv"), row.names = FALSE)
    cat("\nProblem cells saved to: calflora_problem_cells.csv\n")
  } else {
    cat("\n✓ All grid cells downloaded successfully!\n")
  }
  
  # Create visualization
  p <- ggplot(grid_status, aes(x = col, y = row, fill = status)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(
      values = c(
        "ok" = "green3",
        "no_records" = "yellow2",
        "empty_or_corrupt" = "orange",
        "missing" = "red",
        "read_error" = "darkred"
      )
    ) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = "CalFlora Download Grid Status",
      subtitle = sprintf("%d rows × %d cols = %d cells", n_rows, n_cols, n_rows * n_cols),
      x = "Column",
      y = "Row",
      fill = "Status"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )
  
  print(p)
  ggsave("calflora_grid_status.png", p, width = 10, height = 8, dpi = 300)
  cat("\nVisualization saved to: calflora_grid_status.png\n")
  
  return(list(
    grid_status = grid_status,
    problem_cells = problem_cells,
    summary = table(grid_status$status)
  ))
}

# To re-download only the problem cells:
redownload_problem_cells <- function(problem_cells_file = "Temp/calflora_problem_cells.csv",
                                     output_dir = "Data_Raw/Occurrences/calflora_data") {
  
  problem_cells <- read.csv(problem_cells_file)
  
  if (nrow(problem_cells) == 0) {
    cat("No problem cells to redownload!\n")
    return()
  }
  
  # Recreate the geometry
  california_sf <- states(cb = TRUE, resolution = "500k") %>%
    filter(STUSPS == "CA") %>%
    st_make_valid()
  
  california_proj <- st_transform(california_sf, crs = 3310)
  california_buffer <- st_buffer(california_proj, dist = 100000)
  california_buffer_wgs84 <- st_transform(california_buffer, crs = 4326)
  ca_bbox <- st_bbox(california_buffer_wgs84)
  
  total_area <- (ca_bbox["xmax"] - ca_bbox["xmin"]) * (ca_bbox["ymax"] - ca_bbox["ymin"])
  cell_area <- total_area / 500
  cell_side <- sqrt(cell_area)
  
  base_url <- paste0("https://calflora.org/app/download?cols=ID,Plant,Photo,Observer,Date,Source,County,Location+Description,Access,Accuracy:+Square+Meters,Citation,Calrecnum,Common+Name,Cover,Dataset,Date/Time,Datum,Detail+URL,Distribution,Distribution+Code,Elevation+(m),Entered+by,Error+Radius+(m),From,Genus,Geometry+Id,Gross+Area+Count,Gross+Area+Units,GrossArea,Group+Number,Habitat,Herbarium,Index+Date,Infested+Area+Count,Infested+Area+Units,InfestedArea,Large+Photo,Latitude,Location+Quality,Longitude,Management+Status,Modify+Date,Most+Recent,National+Ownership,National+Ownership+Code,Native+Status,Natural+Status,Notes,Number+of+Plants,Observation+Date,Organization,Percent+Cover,Percent+Cover+Midpoint,Phenology,Phenology+Code,Photo+URL,Region,Root,Species,Subspecies,Survey+ID,Taxon&format=CSV&mtk=", mtk_token, "&cell=t&georef=xc&georeferenced=t&ostatus=a&pcount=1&doc_type=rs&natural_status=w&nstatus=2&addnloc=L&cch=t&cnabh=t&wint=r")
  
  cat(sprintf("Re-downloading %d problem cells...\n", nrow(problem_cells)))
  
  for (idx in 1:nrow(problem_cells)) {
    i <- problem_cells$row[idx]
    j <- problem_cells$col[idx]
    
    # Calculate cell bounds
    xmin <- ca_bbox["xmin"] + (j - 1) * cell_side
    xmax <- min(ca_bbox["xmin"] + j * cell_side, ca_bbox["xmax"])
    ymin <- ca_bbox["ymin"] + (i - 1) * cell_side
    ymax <- min(ca_bbox["ymin"] + i * cell_side, ca_bbox["ymax"])
    
    cell_bbox <- c(xmin, ymin, xmax, ymax)
    wkt <- create_wkt_polygon(cell_bbox)
    
    output_file <- file.path(output_dir, sprintf("calflora_row%03d_col%03d.csv", i, j))
    
    cat(sprintf("[%d/%d] Re-downloading cell (%d, %d)...\n", 
                idx, nrow(problem_cells), i, j))
    
    download_calflora_data(wkt, output_file, base_url)
    
    Sys.sleep(2)  # Be extra polite when re-downloading
  }
  
  cat("\nRe-download complete! Run check_calflora_grid() again to verify.\n")
}

#### Download occurrences ####
download_calflora_grid(n_chunks = 500, output_dir = "Data_Raw/Occurrences/calflora_data")

#### Map occurrences to make verify everything correct ####
# Read in all files, merge
calflora_files <- list.files("Data_Raw/Occurrences/calflora_data", full.names = TRUE)
calflora_raw <- map_dfr(calflora_files, ~read_csv(.x, col_types = cols(.default = "c")))

# Convert to numeric
calflora_raw$Latitude  <- suppressWarnings(as.numeric(calflora_raw$Latitude))
calflora_raw$Longitude <- suppressWarnings(as.numeric(calflora_raw$Longitude))

# Set California bounds
xlim <- c(-125, -114)
ylim <- c(32.3, 42.2)

par(mar = c(3, 3, 2, 1))  # smaller margins
plot.new()
plot.window(xlim = xlim, ylim = ylim, asp = 1.3)

# Add points
points(calflora_raw$Longitude, calflora_raw$Latitude,
       pch = 20, cex = 0.25, col = "black")

# Draw red outline *on top*
map("state", "california", xlim = xlim, ylim = ylim,
    fill = FALSE, col = "red", lwd = 1.5, add = TRUE)

# Add concise title
title(main = "Calflora Occurrences", cex.main = 0.9)
box()

#### Check where download may have failed and re-download ####
grid_checked <- check_calflora_grid()

# Run function
redownload_problem_cells()
grid_checked <- check_calflora_grid()

## Map again just to  make sure

# Convert to numeric
calflora_raw$Latitude  <- suppressWarnings(as.numeric(calflora_raw$Latitude))
calflora_raw$Longitude <- suppressWarnings(as.numeric(calflora_raw$Longitude))

# Set California bounds
xlim <- c(-125, -114)
ylim <- c(32.3, 42.2)

# Make plot larger and aspect-correct
par(mar = c(3, 3, 2, 1))  # smaller margins
plot.new()
plot.window(xlim = xlim, ylim = ylim, asp = 1.3)

# Add points
points(calflora_raw$Longitude, calflora_raw$Latitude,
       pch = 20, cex = 0.25, col = "black")

# Draw red outline *on top*
map("state", "california", xlim = xlim, ylim = ylim,
    fill = FALSE, col = "red", lwd = 1.5, add = TRUE)

# Add concise title
title(main = "Calflora Occurrences", cex.main = 0.9)
box()

