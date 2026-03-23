library(data.table)
library(rotl)
library(tidyverse)
library(ape)

#read in the interaction data
interaction_df <- fread("Data_Clean/Interactions/ints_final_clean_slim.csv")

#read in the checklist
checklist_df<- fread("Data_Clean/Species_Checklists/checklist_cleaned.csv")


###After harmonization, we read in the harmonized list and compare it with the species name on the interaction data and the checklist

group_name <- "hoverflies" #change names for other taxon groups

if (group_name == "bees" | group_name == "hoverflies"){
harmonized_names <- fread(paste0("Data_Raw/Phylogeny/", group_name,"_names_tree_harm.csv"))%>%
  select(name_from_tree, names_harm)} else
{harmonized_names<-fread("Data_Raw/Phylogeny/leps_names_tree_harm.csv")%>%
 select(name_from_tree, names_harm)}

group_genera<-checklist_df%>%
  filter(taxon== group_name)%>%
  select(genus)%>%
  unique()%>%
  pull()

checklist_sp<-checklist_df%>%
  filter(taxon==group_name)%>%
  select(genus_species)%>%
  unique()%>%
  arrange()%>%
  pull()

interaction_sp<-interaction_df%>%
  filter(sourceTaxonGenusName %in% group_genera)%>%mutate(genus_species = word(sourceTaxonName_harm, 1,2))%>%
  select(genus_species)%>%unique()%>%arrange()%>%filter(!is.na(genus_species))%>%pull()

####Subset the phylogentic tree from Open Tree of Life####

##create a dataframe with species names and source to input into the tnrs_match_names() function
sp_names_to_input <- data.frame(
  species = c(checklist_sp, interaction_sp), 
  source = c(rep("checklist", length(checklist_sp)), 
             rep("interaction_data", length(interaction_sp)))) %>%
  arrange(species, source) %>%  # This puts "checklist" before "interaction_data" alphabetically
  distinct(species, .keep_all = TRUE)%>%
  left_join(harmonized_names, join_by(species == names_harm))

#double check to make sure any names that overlap between the interaction data and the checklist has the source listed as "checklist"
sp_names_to_input%>%
  filter(species == intersect(interaction_sp, checklist_sp)[1])


#there are species in our data that has multiple matches with names in Open Tree of life, we will priorize the species name that has an exact match (the same name)
duplicated_names_to_keep<-sp_names_to_input%>%
  filter(!is.na(name_from_tree))%>%
  group_by(species)%>%
  filter(n()>1)%>%
  mutate(priority =
           case_when( species == str_replace_all(name_from_tree, "_", " ") ~ 1,
                      TRUE ~ 2)) %>%
  slice_min(priority, n=1, with_ties = FALSE)%>%
  ungroup()%>%
  select(name_from_tree)%>%pull()

#we will remove the other matches
names_to_remove<-sp_names_to_input%>%
  filter(!is.na(name_from_tree))%>%
  group_by(species)%>%
  filter(n()>1)%>%
  filter(!name_from_tree %in% duplicated_names_to_keep)%>%
  ungroup()%>%
  select(name_from_tree)%>%pull()

#final species list to be input into tnrs_match_names()
poll_species <- sp_names_to_input%>%
  filter(!is.na(name_from_tree))%>%
  filter(!name_from_tree %in%  names_to_remove)%>%
  select(name_from_tree)%>%
  pull()

##Now we run tnrs_match_names() to get the species id for extracting the phylogentic tree
#we split into chunks of 1000 names to avoid overloading
chunk_size <- 1000
poll_matches_list <- list()

for(i in seq(1, length(poll_species), by = chunk_size)) {
  chunk <- poll_species[i:min(i + chunk_size - 1, length(poll_species))]
  poll_matches_list[[length(poll_matches_list) + 1]] <- tnrs_match_names(names = chunk)
  Sys.sleep(5)
}

poll_matches <- do.call(rbind, poll_matches_list)

if (group_name == "moths"){
#add this synonym because the prioritzed one cannot be matched with an ott_id
add_on<-tnrs_match_names(names ="Schinia_meskeana")

poll_matches <-rbind(poll_matches, add_on)}

#check to see if any genus has no matching id in the system
sum(is.na(poll_matches$ott_id))

#remove rows that has no match (id = NAs) in the system
good_matches <- poll_matches[!is.na(poll_matches$ott_id),]

#see which row has input name is different from the matched name (search_string !== unique_name)
good_matches%>%
  left_join(sp_names_to_input%>%mutate(name_from_tree = str_to_lower(name_from_tree)), join_by(search_string == name_from_tree))%>%
  mutate(search_string = str_replace_all(search_string, "_", " "),
         unique_name = str_to_lower(unique_name))%>%
  filter(!search_string == unique_name)

#count the numbers of unique id that are duplicated after filtered by search_string == unique_name
good_matches%>%
  left_join(sp_names_to_input%>%mutate(name_from_tree = str_to_lower(name_from_tree)), join_by(search_string == name_from_tree))%>%
  mutate(search_string = str_replace_all(search_string, "_", " "),
         unique_name = str_to_lower(unique_name))%>%
  filter(search_string == unique_name)%>%
  group_by(ott_id) %>%
  filter(n() > 1) 


#If any rows have search_string !== unique_name and there are duplicated ids, we want to filter out the rows having search_string !== unique_name
final_id_to_input<-good_matches%>%
  left_join(sp_names_to_input%>%mutate(name_from_tree = str_to_lower(name_from_tree)), join_by(search_string == name_from_tree))%>%
  mutate(search_string = str_replace_all(search_string, "_", " "),
         unique_name = str_to_lower(unique_name))%>%
  filter(search_string == unique_name)

#make sure there is no duplicated id 
sum(duplicated(final_id_to_input$ott_id))

##Check which OTT IDs are actually in the tree
# This step filters out problematic IDs
ott_ids_to_check <- final_id_to_input$ott_id

# Test IDs in small batches to find problematic ones
valid_id_df <- data.frame() #tracking the valid IDs
invalid_id_df<-data.frame() #tracking the problematic IDs

for(i in 1:length(ott_ids_to_check)) {
  
  test_id <- ott_ids_to_check[i]
  # Try to get a tree with just this ID and one known good ID
  # 3272412 is halictus lucidipennis  - a reliable reference
  try_tree <- try(tol_induced_subtree(ott_ids = c(3272412, test_id)), silent = TRUE)
  
  if(!inherits(try_tree, "try-error")) {

    valid_df<-data.frame(species = final_id_to_input$unique_name[final_id_to_input$ott_id == test_id], 
      id= test_id, 
      name_from_tree = final_id_to_input$search_string[final_id_to_input$ott_id == test_id])

    valid_id_df <- rbind(valid_id_df, valid_df)
  } else {

    cat("Saving problematic ID:", test_id, "-", final_id_to_input$unique_name[final_id_to_input$ott_id == test_id], "\n")
    
    invalid_df<-data.frame(species = final_id_to_input$unique_name[final_id_to_input$ott_id == test_id], 
      id= test_id, 
      name_from_tree = final_id_to_input$search_string[final_id_to_input$ott_id == test_id])
    
    invalid_id_df <- rbind(invalid_id_df, invalid_df)
    
  }
}

#we want to make sure the species with invalid ids do not have another synonym that actually has an id
if (nrow(invalid_id_df) > 0){
names_to_search_again<-invalid_id_df%>%
  left_join(sp_names_to_input%>%mutate(name_from_tree = str_to_lower(str_replace_all(name_from_tree, "_", " "))), join_by(name_from_tree == name_from_tree))%>%
  select(species.y)%>%
  left_join(harmonized_names, join_by(species.y == names_harm))

names_to_add<-tnrs_match_names(names = names_to_search_again$name_from_tree)

#filter out the invalid id
names_to_add_df<-names_to_add%>%
  filter(!ott_id %in% invalid_id_df$id)

ott_ids_to_check <- c(valid_id_df$id, names_to_add_df$ott_id)

for(i in 1:length(ott_ids_to_check)) {
  test_id <- ott_ids_to_check[i]
  try_tree <- try(tol_induced_subtree(ott_ids = c(3272412, test_id)), silent = TRUE)
  if(!inherits(try_tree, "try-error")) {
    valid_df<-data.frame(species = names_to_add$unique_name[names_to_add$ott_id == test_id], id= test_id, name_from_tree = names_to_add$search_string[names_to_add$ott_id == test_id])
    valid_id_df <- rbind(valid_id_df, valid_df)
  } else {
    cat("Saving problematic ID:", test_id, "-",
        names_to_add$unique_name[names_to_add$ott_id == test_id], "\n")
    invalid_df<-data.frame(species = names_to_add$unique_name[names_to_add$ott_id == test_id], id= test_id, name_from_tree = names_to_add$search_string[names_to_add$ott_id == test_id])
    invalid_id_df <- rbind(invalid_id_df, invalid_df)
    
  }
}
}

### Now we build the tree with only valid IDs
poll_tree <- tol_induced_subtree(ott_ids = valid_id_df$id, label_format ="id")

data_to_match<-rbind(final_id_to_input, names_to_add_df%>%
  left_join(sp_names_to_input%>%mutate(name_from_tree= str_to_lower(name_from_tree)), join_by(search_string == name_from_tree))%>%
  mutate(search_string =  str_replace_all(search_string, "_", " ")))

# Now we match the tip labels with the name_from_tree using the ott_id
for(i in 1:length(poll_tree$tip.label)) {
  
  # Extract OTT ID from tip label
  ott_id_str <- sub(".*ott", "", poll_tree$tip.label[i])  
  
  if(ott_id_str %in% data_to_match$ott_id) {
    
    # Match it with the original name
    original_name <- data_to_match$search_string[data_to_match$ott_id == ott_id_str]
    
    poll_tree$tip.label[i]<-original_name
    
    print(original_name)
  }
}

#we need to drop and select only one tips label if any are pointing to the same species after the harmonization
duplicated_names_to_keep<-data.frame(tip_label= poll_tree$tip.label)%>%
       left_join(data_to_match, join_by(tip_label == search_string))%>%
  mutate(species = str_to_lower(species))%>%
       group_by(species)%>%
       filter(n()>1)%>%
      mutate(priority =
               case_when( species == tip_label ~ 1,
                          score == max(score) & species != tip_label ~2, 
                          TRUE ~ 3)) %>%
  slice_min(priority, n=1, with_ties = FALSE)%>%
  ungroup()%>%
  select(tip_label)%>%pull()
    

#Drop the tips that is dupulicated
tips_to_drop <- data.frame(tip_label= poll_tree$tip.label)%>%
  left_join(final_id_to_input, join_by( tip_label == search_string))%>%
  mutate(species = str_to_lower(species))%>%
  group_by(species)%>%
  filter(n()>1)%>%
  ungroup()%>%
  filter(! tip_label %in% duplicated_names_to_keep)%>%
  select(tip_label)%>%pull()

pruned_tree <- drop.tip(poll_tree, tips_to_drop)

#Add branch lengths
final_tree <- compute.brlen(pruned_tree, method = "Grafen", power = 0.5)

#Update the tip labels with the harmonized names
for(i in 1:length(final_tree$tip.label)) {
  
  #look for the harmonized name
  search_name<-final_tree$tip.label[i]
  
  if(search_name %in% final_id_to_input$search_string) {
    
    # Match it with the original name
    original_name <- final_id_to_input$species[final_id_to_input$search_string == search_name]
    
    final_tree$tip.label[i]<-original_name
    
    print(original_name)
  }
}

#check to see if the length of input is the same as the output
length(final_tree$tip.label) == (length(poll_tree$tip.label)-length(tips_to_drop))

#save the tree
write.tree(final_tree, paste0("Data_Raw/Phylogeny/", group_name, "_tree.nwk"))


## Create a summary table to summarize species that were not on the tree

#names that are only at genus level is not included
missing_genera<-interaction_df%>%
  filter(sourceTaxonGenusName %in% group_genera)%>%
  mutate(genus_species = word(sourceTaxonName_harm, 1,2))%>%
  filter(is.na(genus_species))%>%
  select(sourceTaxonName_harm)%>%
  unique()%>%
  rename(species = sourceTaxonName_harm)%>%
  mutate(source = "interaction_data",
         reason = "name at genus level")

#names that has no match with the harmonized taxonomy
sp_with_no_match<-sp_names_to_input%>%
  filter(is.na(name_from_tree))%>%
  select(species, source)%>%
  mutate(reason = "no match in the harmonized taxonomy")
  

#names that are not avaliable on the tree
names_not_on_the_tree<-invalid_id_df%>%
  select(name_from_tree)%>%
  left_join(sp_names_to_input%>%mutate(name_from_tree = str_to_lower(str_replace_all(name_from_tree, "_", " "))))%>%
  select(species, source)%>%
  mutate(reason = "names not avaliable on the tree")

  
#the summary dataframe
missing_names_summary<-rbind(missing_genera, sp_with_no_match, names_not_on_the_tree)

#Save the summary
write.csv(missing_names_summary, paste0("Data_Raw/Phylogeny/", group_name, "_without_phylogeny.csv"), row.names = FALSE)