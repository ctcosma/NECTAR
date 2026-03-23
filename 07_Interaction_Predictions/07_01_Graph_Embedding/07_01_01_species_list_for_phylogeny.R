#### Load Packages ####
library(data.table)
library(rotl)
library(tidyverse)
library(ape)

#the interaction data
interaction_df <- fread("Data_Clean/Interactions/ints_final_clean_slim.csv")

#master checklist (all native species)
checklist_df<- fread("Data_Clean/Species_Checklists/checklist_cleaned.csv")

###We need to harmonize the species names in Open Tree of Life to maximize the number of species matches between our checklist and their phylogentic data

##To do so we need to extract all species names at a higher taxonmoic level for a given taxon group from the Open Tree Taxonomy
#leps
lep_match <- tnrs_match_names("Lepidoptera")
leps_tree<-taxonomy_subtree(ott_id = lep_match$ott_id, output_format = "phylo", label_format="name")

leps_names_from_the_tree<-data.frame(name_from_tree = leps_tree$tip.label)%>%
  mutate(species =  str_replace_all(name_from_tree, "_", " "))%>%
  mutate(species = str_to_sentence(species))%>%
  mutate(sp_name = word(species, 2))%>%
  filter(!is.na(sp_name))%>%
  mutate(genus_species = word(species, 1, 2))%>%
  filter(!grepl("\\(genus", name_from_tree))%>%
  filter(!grepl("\\(inconsistent", name_from_tree))%>%
  filter(!grepl("\\(merged", name_from_tree))%>%
  filter(!grepl("unclassified", name_from_tree))%>%
  mutate(genus_species = str_remove_all(genus_species, "[\\(\\)']"))%>%
  select(name_from_tree, genus_species)

write.csv(leps_names_from_the_tree, "Data_Raw/Phylogeny/leps_names_from_raw_tree.csv", row.names = FALSE)

#bees
bees_match <- tnrs_match_names("Apoidea")
bees_tree<-taxonomy_subtree(ott_id = bees_match$ott_id, output_format = "phylo", label_format="name")

bees_names_from_the_tree<-data.frame(name_from_tree = bees_tree$tip.label)%>%
  mutate(species =  str_replace_all(name_from_tree, "_", " "))%>%
  mutate(species = str_to_sentence(species))%>%
  mutate(sp_name = word(species, 2))%>%
  filter(!is.na(sp_name))%>%
  mutate(genus_species = word(species, 1, 2))%>%
  filter(!grepl("\\(genus", name_from_tree))%>%
  filter(!grepl("\\(inconsistent", name_from_tree))%>%
  filter(!grepl("\\(merged", name_from_tree))%>%
  filter(!grepl("unclassified", name_from_tree))%>%
  mutate(genus_species = str_remove_all(genus_species, "[\\(\\)']"))%>%
  select(name_from_tree, genus_species)

write.csv(bees_names_from_the_tree, "Data_Raw/Phylogeny/bees_names_from_raw_tree.csv", row.names = FALSE)

#hoverflies
hoverflies_match <- tnrs_match_names("Syrphidae")
hoverflies_tree<-taxonomy_subtree(ott_id = hoverflies_match$ott_id, output_format = "phylo", label_format="name")

hoverflies_names_from_the_tree<-data.frame(name_from_tree = hoverflies_tree$tip.label)%>%
  mutate(species =  str_replace_all(name_from_tree, "_", " "))%>%
  mutate(species = str_to_sentence(species))%>%
  mutate(sp_name = word(species, 2))%>%
  filter(!is.na(sp_name))%>%
  mutate(genus_species = word(species, 1, 2))%>%
  filter(!grepl("\\(genus", name_from_tree))%>%
  filter(!grepl("\\(inconsistent", name_from_tree))%>%
  filter(!grepl("\\(merged", name_from_tree))%>%
  filter(!grepl("unclassified", name_from_tree))%>%
  mutate(genus_species = str_remove_all(genus_species, "[\\(\\)']"))%>%
  select(name_from_tree, genus_species)

write.csv(hoverflies_names_from_the_tree, "Data_Raw/Phylogeny/hoverflies_names_from_raw_tree.csv", row.names = FALSE)


##### Harmonize names of open tree #####
source("Code/99_Supporting/name_harmonization.R")

# Read in names from opentree
bee_names_tree <- read.csv("Data_Raw/Phylogeny/bees_names_from_raw_tree.csv")
lep_names_tree <- read.csv("Data_Raw/Phylogeny/leps_names_from_raw_tree.csv")
hoverflies_names_tree <- read.csv("Data_Raw/Phylogeny/hoverflies_names_from_raw_tree.csv")

names_bees <- bee_names_tree %>% pull(genus_species) %>% unique()
names_leps <- lep_names_tree %>% pull(genus_species) %>% unique()
names_hoverflies <- hoverflies_names_tree %>% pull(genus_species) %>% unique()

# Harmonize names
names_harm_bees <- harmonize_names(
  names = names_bees,
  higher_tax = "Hymenoptera",
  names_file = "Temp/names_opentree_bees.tsv",
  names_harm_file = "Temp/names_harmonized_opentree_bees.tsv"
)

names_harm_leps <- harmonize_names(
  names = names_leps,
  higher_tax = "Lepidoptera",
  names_file = "Temp/names_opentree_leps.tsv",
  names_harm_file = "Temp/names_harmonized_opentree_leps.tsv"
)

names_harm_hoverflies <- harmonize_names(
  names = names_hoverflies,
  higher_tax = "Diptera",
  names_file = "Temp/names_opentree_hoverflies.tsv",
  names_harm_file = "Temp/names_harmonized_opentree_hoverflies.tsv"
)

# Append harmonized names back to original tree name datasets
bee_names_tree_harm <- bee_names_tree %>%
  left_join(names_harm_bees, by = join_by("genus_species" == "ScientificName"))
lep_names_tree_harm <- lep_names_tree %>%
  left_join(names_harm_leps, by = join_by("genus_species" == "ScientificName"))
hoverflies_names_tree_harm <- hoverflies_names_tree %>%
  left_join(names_harm_hoverflies, by = join_by("genus_species" == "ScientificName"))

# Get higher taxonomy
bee_names_tree_harm <- update_taxonomy(bee_names_tree_harm, tax = "animal")
lep_names_tree_harm <- update_taxonomy(lep_names_tree_harm, tax = "animal")
hoverflies_names_tree_harm <- update_taxonomy(hoverflies_names_tree_harm, tax = "animal")

# Save
write.csv(bee_names_tree_harm, "Data_Raw/Phylogeny/bees_names_tree_harm.csv")
write.csv(lep_names_tree_harm, "Data_Raw/Phylogeny/leps_names_tree_harm.csv")
write.csv(hoverflies_names_tree_harm, "Data_Raw/Phylogeny/hoverflies_names_tree_harm.csv")