#### Load packages ####
library(data.table)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)

#### Read Data ####
#the interaction data
interaction_df<-fread("Data_Clean/Interactions/ints_final_clean_withCrops_slim.csv")%>%
  filter(is.na(isCrop))%>%
  mutate(source_genus_species = word(sourceTaxonName_harm,1,2))

checklist_df<- fread("Data_Clean/Species_Checklists/checklist_cleaned.csv")

#Var_explained= 1.0
#read in the interaction network of the checklist
bee_net_1.0<-fread("Temp/bees_var_exp_1.0_prediction_summary.csv")

hoverflies_net_1.0<-fread("Temp/hoverflies_var_exp_1.0_prediction_summary.csv")

butterflies_net_1.0<-fread("Temp/butterflies_var_exp_1.0_prediction_summary.csv")

moths_net_1.0<-fread("Temp/moths_var_exp_1.0_prediction_summary.csv")

#count the number of observed interaction in each taxon group
bee_total_obs_interaction_1.0<-bee_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1))%>%
  summarise(total_obs_interaction = sum(observed_interaction))%>%pull()

hoverflies_total_obs_interaction_1.0<-hoverflies_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1))%>%
  summarise(total_obs_interaction = sum(observed_interaction))%>%pull()

butterflies_total_obs_interaction_1.0<-butterflies_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1))%>%
  summarise(total_obs_interaction = sum(observed_interaction))%>%pull()

moths_total_obs_interaction_1.0<-moths_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1))%>%
  summarise(total_obs_interaction = sum(observed_interaction))%>%pull()

#Count the number of interaction gained in each type of learning after prediction
bee_total_gain_1.0<-bee_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"))%>%
  group_by(type_of_learning)%>%
  filter(predicted_interaction > observed_interaction)%>%
  summarise(total_gain = sum(predicted_interaction))%>%ungroup()%>%
  mutate(group = "bees")%>%
  mutate(total_sp_used = length(unique(bee_net_1.0$genus_species)),
         total_gain_per_sp = total_gain/total_sp_used,
         total_checklist_sp = checklist_df%>%filter(taxon =="bees")%>%
           select(genus_species)%>%
           pull()%>%unique()%>%length(),
         phylogenetic_coverage = total_sp_used/total_checklist_sp)%>%
  select(group, type_of_learning, total_checklist_sp, phylogenetic_coverage, total_sp_used, total_gain, total_gain_per_sp)

hoverflies_total_gain_1.0<-hoverflies_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"))%>%
  group_by(type_of_learning)%>%
  filter(predicted_interaction > observed_interaction)%>%
  summarise(total_gain = sum(predicted_interaction))%>%ungroup()%>%
  mutate(group = "hoverflies")%>%
  mutate(total_sp_used = length(unique(hoverflies_net_1.0$genus_species)),
         total_gain_per_sp = total_gain/total_sp_used,
         total_checklist_sp = checklist_df%>%filter(taxon =="hoverflies")%>%
           select(genus_species)%>%
           pull()%>%unique()%>%length(),
         phylogenetic_coverage = total_sp_used/total_checklist_sp)%>%
  select(group, type_of_learning, total_checklist_sp, phylogenetic_coverage, total_sp_used, total_gain, total_gain_per_sp)

butterflies_total_gain_1.0<-butterflies_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"))%>%
  group_by(type_of_learning)%>%
  filter(predicted_interaction > observed_interaction)%>%
  summarise(total_gain = sum(predicted_interaction))%>%ungroup()%>%
  mutate(group = "butterflies")%>%
  mutate(total_sp_used = length(unique(butterflies_net_1.0$genus_species)),
         total_gain_per_sp = total_gain/total_sp_used,
         total_checklist_sp = checklist_df%>%filter(taxon =="butterflies")%>%
           select(genus_species)%>%
           pull()%>%unique()%>%length(),
         phylogenetic_coverage = total_sp_used/total_checklist_sp)%>%
  select(group, type_of_learning, total_checklist_sp, phylogenetic_coverage, total_sp_used, total_gain, total_gain_per_sp)

moths_total_gain_1.0<-moths_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"))%>%
  group_by(type_of_learning)%>%
  filter(predicted_interaction > observed_interaction)%>%
  summarise(total_gain = sum(predicted_interaction))%>%ungroup()%>%
  mutate(group = "moths")%>%
  mutate(total_sp_used = length(unique(moths_net_1.0$genus_species)),
         total_gain_per_sp = total_gain/total_sp_used,
         total_checklist_sp = checklist_df%>%filter(taxon =="moths")%>%
           select(genus_species)%>%
           pull()%>%unique()%>%length(),
         phylogenetic_coverage = total_sp_used/total_checklist_sp)%>%
  select(group, type_of_learning, total_checklist_sp, phylogenetic_coverage, total_sp_used, total_gain, total_gain_per_sp)



#look at the genus level interaction gain
bee_interaction_gain_1.0<-bee_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"),
    pollinator_genus = word(genus_species , 1))%>%
  group_by(pollinator_genus)%>%
  mutate(n_sp_in_genus = n_distinct(genus_species))%>%
  ungroup()%>%
  filter(predicted_interaction > observed_interaction)%>%
  group_by(pollinator_genus, type_of_learning)%>%
  summarise(total_interaction_gained = sum(predicted_interaction),
            n_sp_in_genus = first(n_sp_in_genus))%>%
  ungroup()%>%
  mutate(interaction_gained_per_sp_genus = total_interaction_gained / n_sp_in_genus)

hoverflies_interaction_gain_1.0<-hoverflies_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"),
    pollinator_genus = word(genus_species , 1))%>%
  group_by(pollinator_genus)%>%
  mutate(n_sp_in_genus = n_distinct(genus_species))%>%
  ungroup()%>%
  filter(predicted_interaction > observed_interaction)%>%
  group_by(pollinator_genus, type_of_learning)%>%
  summarise(total_interaction_gained = sum(predicted_interaction),
            n_sp_in_genus = first(n_sp_in_genus))%>%
  ungroup()%>%
  mutate(interaction_gained_per_sp_genus = total_interaction_gained / n_sp_in_genus)


butterflies_interaction_gain_1.0<-butterflies_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"),
    pollinator_genus = word(genus_species , 1))%>%
  group_by(pollinator_genus)%>%
  mutate(n_sp_in_genus = n_distinct(genus_species))%>%
  ungroup()%>%
  filter(predicted_interaction > observed_interaction)%>%
  group_by(pollinator_genus, type_of_learning)%>%
  summarise(total_interaction_gained = sum(predicted_interaction),
            n_sp_in_genus = first(n_sp_in_genus))%>%
  ungroup()%>%
  mutate(interaction_gained_per_sp_genus = total_interaction_gained / n_sp_in_genus)


moths_interaction_gain_1.0<-moths_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"),
    pollinator_genus = word(genus_species , 1))%>%
  group_by(pollinator_genus)%>%
  mutate(n_sp_in_genus = n_distinct(genus_species))%>%
  ungroup()%>%
  filter(predicted_interaction > observed_interaction)%>%
  group_by(pollinator_genus, type_of_learning)%>%
  summarise(total_interaction_gained = sum(predicted_interaction),
            n_sp_in_genus = first(n_sp_in_genus))%>%
  ungroup()%>%
  mutate(interaction_gained_per_sp_genus = total_interaction_gained / n_sp_in_genus)

#the plots
# summary_plot_1.0<- rbind(bee_total_gain_1.0, hoverflies_total_gain_1.0, butterflies_total_gain_1.0, moths_total_gain_1.0)%>%
#   rename(type = type_of_learning, total = total_gain)%>%
#   rbind(data.frame(group = c("bees", "hoverflies", "butterflies", "moths"), type = c(rep("observed interaction", 5)), total= c(bee_total_obs_interaction_1.0, hoverflies_total_obs_interaction_1.0, butterflies_total_obs_interaction_1.0, moths_total_obs_interaction_1.0)))%>%
#   ggplot(aes(x = group, y = total, fill = type)) +
#   geom_col(position = "stack") +
#   coord_flip()+
#   theme_cowplot()+
#   theme(legend.position = "bottom",
#         axis.text.y = element_text(size=9),
#         axis.title.x =  element_text(size=9))+
#   labs(x= "Group", y = "No. of current interaction", fill = "Type")
# 
# ggsave2(paste0("~/Documents/research/Morpho_Interactions/Code/07_Interaction_Prediction/07_05_graph_embedding/output/plots/interaction_gain_summary_1.0.png"), 
#         summary_plot_1.0 , dpi = 300, bg = "white", width = 8, height = 5)
# 

#output a summary table
summary_df<-rbind(bee_total_gain_1.0, hoverflies_total_gain_1.0, butterflies_total_gain_1.0, moths_total_gain_1.0)
write.csv(summary_df, "Data_Clean/Graph_Embedding/predicted_summary.csv", row.names = FALSE)


#save the actual interaction prediction in cvs format for downstream analyses
bees_net<-bee_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"))

hoverflies_net<-hoverflies_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"))

butterflies_net<-butterflies_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"))

moths_net<-moths_net_1.0%>%
  mutate(observed_interaction= case_when(
    observed_interaction == "FALSE" ~ 0,
    TRUE ~ 1),
    predicted_interaction= case_when(
      predicted_interaction == "FALSE" ~0,
      TRUE ~ 1),
    type_of_learning = case_when(
      genus_species %in% interaction_df$source_genus_species ~ "In-sample learning",
      TRUE ~ "Out-of-sample learning"))

write.csv(bees_net, "Data_Clean/Graph_Embedding/Results_Long/bees_predicted_net.csv", row.names = FALSE)

write.csv(hoverflies_net, "Data_Clean/Graph_Embedding/Results_Long/hoverflies_predicted_net.csv", row.names = FALSE)

write.csv(butterflies_net, "Data_Clean/Graph_Embedding/Results_Long/butterflies_predicted_net.csv", row.names = FALSE)

write.csv(moths_net, "Data_Clean/Graph_Embedding/Results_Long/moths_predicted_net.csv", row.names = FALSE)


