#### Load Packages ####
library(bdc)
library(tidyverse)

#### Prepare BDC Metadata file ####
# metadata file is a csv spreadsheet, first two columns have to be "datasetName" 
# (e.g., "gbif_formatted") and "fileName" (e.g., "gbif_formatted.csv"), then each 
# other column is the column name you want (darwin core format), and its name in 
# each corresponding dataset, and they will all be harmonized to the column name

# Load subsets of the occurrence data to prepare bdc metadata file
gbif = read.csv("Data_Raw/Formatted/gbif_formatted.csv", nrows = 100)
ecdysis = read.csv("Data_Raw/Formatted/ecdysis_formatted.csv", nrows = 100)
chesshire = read.csv("Data_Raw/Formatted/chesshire_formatted.csv", nrows = 100)
calflora = read.csv("Data_Raw/Formatted/calflora_formatted.csv", nrows = 100)
cch2 = read.csv("Data_Raw/Formatted/cch2_formatted.csv", nrows = 100)
cch2_pheno = read.csv("Data_Raw/Formatted/cch2_pheno_formatted.csv", nrows = 100)
mpg = read.csv("Data_Raw/Formatted/mpg_formatted.csv", nrows = 100)
idigbio = read.csv("Data_Raw/Formatted/idigbio_formatted.csv", nrows = 100)


# gbif, ecdysis, chesshire, cch2 all in Darwincore
# Find matching column names between them
datasets <- list(gbif, ecdysis, cch2, chesshire, cch2_pheno)
matching_columns <- Reduce(intersect, lapply(datasets, names))
matching_columns <- data.frame(t(matching_columns))

# save and manually prepare the bdc_meta.csv file in spreadsheet program
write.csv(matching_columns, "Data_Raw/Formatted/matching_columns.csv", row.names = F)

# Be sure to add Phenology and Phenology.Code for calflora data

# Calflora, mpg, and idigbio are the weird ones, have to manually add the corresponding columns in
names(mpg)
names(calflora)
names(idigbio)

##### Dataset Harmonization #####

# Clear workspace
rm(list = ls(all = TRUE))
gc()

# Load the required metadata file
bdc_meta <- read_csv("Data_Raw/Formatted/bdc_meta.csv")

# Merge and standardize the datasets with bdc
database <-
  bdc_standardize_datasets(
    metadata = bdc_meta,
    format = "csv",
    overwrite = TRUE,
    save_database = FALSE) %>%
  #remove non UTF-8 characters
  mutate(scientificName = iconv(scientificName, from = "latin1", to = "UTF-8", sub = ""))

#### Save output harmonized database ####
write.csv(database, "Data_Raw/Formatted/database.csv", row.names = F)
