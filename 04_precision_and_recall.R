library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(forcats) 
library(ggpubr)
library(glue)
library(reshape2)
library(phyloseq)

text_size = 7
custom_theme = function(){ 
  theme_classic() %+replace% 
    theme(
      #text elements
      plot.title = element_text(size = text_size),                #set font size
      plot.subtitle = element_text(size = text_size),               #font size
      plot.caption = element_text(size = text_size),               #right align
      axis.title = element_text(size = text_size),               #font size
      axis.text = element_text(size = text_size),                #font size
      axis.text.x = element_text(size = text_size),
      axis.text.y = element_text(size = text_size), 
      legend.text = element_text(size = text_size), 
      legend.title = element_text(size = text_size, hjust=0), 
      strip.text = element_text(size = text_size),
      strip.background = element_blank(), 
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
      panel.grid.minor = element_blank(), 
      #panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid', colour = "grey")
      panel.grid.major = element_blank()
    )
}

##############################################################################
############################## Load csv files ################################  
##############################################################################

raw12s = read.csv("/outputs/12s_allmethods.csv")
rawcoi = read.csv("/outputs/coi_allmethods.csv")

raw12s$Marker_gene = "12S"
rawcoi$Marker_gene = "COI"

df01 = rbind(raw12s, rawcoi)
df01$Method = paste(df01$Marker_gene, df01$Method, sep = " - ")

df01$Sample = gsub("12s-|coi-", "", df01$Sample)

confspecies = c("Specie bootstrap < 70", "Percentage ID < 97%")
df01$Confidence_categories = ifelse(df01$Species %in% confspecies, 
                                    "Confidence score below threhold", 
                                    "Confidence score above threhold")

df01$Sample = gsub("^even", "Even distribution of reads", df01$Sample)
df01$Sample = gsub("staggered", "Staggered distribution of reads", df01$Sample)
df01$Sample = gsub("allalberta", "All Alberta fish species", df01$Sample)
df01$Sample = gsub("blindspot", "Fish identified as taxonomic blindspot", df01$Sample)
df01$Sample = gsub("lowrelat", "Community of species with low relatedness", df01$Sample)
df01$Sample = gsub("mix", "All Alberta and other species of fish", df01$Sample)
df01$Sample = gsub("richplusevenplus", "Blindspot - High richness and high evenness", df01$Sample)
df01$Sample = gsub("richplusevenlow", "Blindspot - High richness and low evenness", df01$Sample)
df01$Sample = gsub("richlowevenlow", "Blindspot - Low richness and low evenness", df01$Sample)
df01$Sample = gsub("richlowevenplus", "Blindspot - Low richness and high evenness", df01$Sample)
df01$Species = gsub("Hemichromis bimaculatus", "Rubricatochromis bimaculatus", df01$Species)

##############################################################################
############# Calculate precision and recall for each samples  ###############  
##############################################################################

############################ EMPIRICAL DATA ##################################

# Load list of Alberta fresh water fishes for common name  
alberta_fish = read.csv("/database/fishbase_alberta_with_missing.csv", header = TRUE, check.names = FALSE)
# Create list to ID fish not from Alberta
list_alberta_fish = unique(alberta_fish$Species)
extra_species = c("Oncorhynchus mykiss x Salmo salar", "Cottus bairdii:cognatus", 
                  "Salvelinus alpinus:malma", "Coregonus artedi:zenithicus")
list_alberta_species = c(list_alberta_fish, extra_species)

df_empraw = subset(df01, grepl("00", df01$Sample) &
                     !grepl("<", df01$Species))

df_empraw$Prediction_classification = ifelse(!df_empraw$Species %in% list_alberta_species, "False positive", "True positive")

df_grouped = df_empraw %>% 
  group_by(Method, Sample, Prediction_classification) %>%
  summarise(Count = n())

df_grouped$GroupingID = paste(df_grouped$Method, df_grouped$Sample, sep = "_")

df_unmelt = dcast(df_grouped, GroupingID ~ Prediction_classification, value.var = "Count", fill = 0)

df_unmelt$Precision = df_unmelt$`True positive` /
  (df_unmelt$`True positive` + df_unmelt$`False positive`)

df_empirical = df_unmelt %>% separate(GroupingID, sep = "_", c("Method", "Sample"))

# Proportion classified as Alberta fish, Non-fish, Non-Alberta fish

df_empraw$Group_taxonomy = ifelse(grepl("Non-", df_empraw$Species), df_empraw$Species, "Alberta fish")

df_empsum_group = df_empraw %>%
  group_by(Method, Sample, Group_taxonomy) %>%
  summarise(Total_asv = sum(Abundance))

df_empproportion = df_empsum_group %>%
  group_by(Method, Sample) %>%
  mutate(Proportion = Total_asv / sum(Total_asv) * 100)

############################ BLINDSPOT COMMUNITY ##################################

# Load list of species included in mock communities
trueblind = read.table("/mock_communities/names_blindpost.txt", sep = ";")
names(trueblind) = c("ID","CellularOrgamism", "Superkingdom", "Kingdom", "Phylum", 
                     "Class", "Order", "Family", "Genus","Species")
trueblind$Species = gsub("_", " ", trueblind$Species)
trueblind$Prediction_classification = "True positive"
trueblind = select(trueblind, c("Species", "Prediction_classification"))
trueblind = trueblind[!duplicated(trueblind), ]

# Extract sample to evaluate 
df_blindsplot01 = subset(df01, Sample == "Fish identified as taxonomic blindspot" & 
                           !grepl("<", df01$Species))

df_blindsplot02 = df_blindsplot01 %>% 
  group_by(Method, Species) %>%
  summarize(Total = n())

listdf = list()
i = 0 

for(method in unique(df_blindsplot02$Method)){
  i = i + 1 
  sub = subset(df_blindsplot02, Method == method)
  combdf = merge(trueblind, sub, by = "Species", all = TRUE)
  combdf = combdf %>%  fill(Method)
  listdf[[i]] = combdf
}

df_blindsplot03 = do.call("rbind", listdf)
df_blindsplot03$Prediction_classification = df_blindsplot03$Prediction_classification %>% 
  replace_na("False positive")
df_blindsplot03$Prediction_classification = ifelse(is.na(df_blindsplot03$Total),
                                                   "False negative",
                                                   df_blindsplot03$Prediction_classification)

df_blindsplot04 = df_blindsplot03 %>% 
  group_by(Method, Prediction_classification) %>%
  summarise(Count = n())

df_blindsplot05 = dcast(df_blindsplot04, Method ~ Prediction_classification, value.var = "Count", fill = 0)

df_blindsplot05$Precision = df_blindsplot05$`True positive` /
  (df_blindsplot05$`True positive` + df_blindsplot05$`False positive`)

df_blindsplot05$Recall = df_blindsplot05$`True positive` /
  (df_blindsplot05$`True positive` + df_blindsplot05$`False negative`)

df_blindsplot05$Sample = "Fish identified as taxonomic blindspot"
df_blindsplot03$Sample = "Fish identified as taxonomic blindspot"

############################ LOW RELATEDNESS ##################################

# Load list of species included in mock communities
truelow = read.table("/mock_communities/names_lowrelat.txt", sep = ";")
names(truelow) = c("ID","CellularOrgamism", "Superkingdom", "Kingdom", "Phylum", 
                   "Class", "Order", "Family", "Genus","Species")
truelow$Species = gsub("_", " ", truelow$Species)
truelow$Prediction_classification = "True positive"
truelow = select(truelow, c("Species", "Prediction_classification"))
truelow = truelow[!duplicated(truelow), ]

# Extract sample to evaluate 
df_low01 = subset(df01, Sample == "Community of species with low relatedness" & 
                    !grepl("<", df01$Species))

df_low02 = df_low01 %>% 
  group_by(Method, Species) %>%
  summarize(Total = n())

listdf = list()
i = 0 

for(method in unique(df_low02$Method)){
  i = i + 1 
  sub = subset(df_low02, Method == method)
  combdf = merge(truelow, sub, by = "Species", all = TRUE)
  combdf = combdf %>%  fill(Method)
  listdf[[i]] = combdf
}

df_low03 = do.call("rbind", listdf)

df_low03$Prediction_classification = df_low03$Prediction_classification %>% 
  replace_na("False positive")
df_low03$Prediction_classification = ifelse(is.na(df_low03$Total),
                                            "False negative",
                                            df_low03$Prediction_classification)
df_low04 = df_low03 %>% 
  group_by(Method, Prediction_classification) %>%
  summarise(Count = n())

df_low05 = dcast(df_low04, Method ~ Prediction_classification, value.var = "Count", fill = 0)

df_low05$Precision = df_low05$`True positive` /
  (df_low05$`True positive` + df_low05$`False positive`)

df_low05$Recall = df_low05$`True positive` /
  (df_low05$`True positive` + df_low05$`False negative`)

df_low05$Sample = "Community of species with low relatedness"
df_low03$Sample = "Community of species with low relatedness"

############################ ALL ALBERTA ##################################

# Because some species didn't have a sequence for each marker gene, they are processed separately. 

# List of marker genes 
list_marker_gene = c("12s", "coi")
list_truealberta = list.files("/mock_communities", 
                              pattern = "albertaspecies.txt", full.names = TRUE)
list_truedf = list()
i = 0 

for(genes in list_marker_gene){
  i = i + 1 
  marker_gene = list_marker_gene[i]
  true = read.table(file = grep(marker_gene, list_truealberta, value = TRUE), sep = ";")
  names(true) = c("ID","CellularOrgamism", "Superkingdom", "Kingdom", "Phylum", 
                  "Class", "Order", "Family", "Genus","Species")
  true$Species = gsub("_", " ", true$Species)
  true$Prediction_classification = "True positive"
  true = select(true, c("Species", "Prediction_classification"))
  list_truedf[[i]] = true
  names(list_truedf)[[i]] = marker_gene
} 

# Extract sample to evaluate 
df_allalberta01 = subset(df01, Sample == "All Alberta fish species" & 
                           !grepl("<", df01$Species))

df_allalberta02 = df_allalberta01 %>% 
  group_by(Method, Species) %>%
  summarize(Total = n())

# Process 12S 

df_allalberta12s = subset(df_allalberta02, grepl("12S", df_allalberta02$Method))

listdf = list()
i = 0 
for(method in unique(df_allalberta12s$Method)){
  i = i + 1 
  sub = subset(df_allalberta12s, Method == method)
  combdf = merge(list_truedf$`12s`, sub, by = "Species", all = TRUE)
  combdf = combdf %>%  fill(Method)
  listdf[[i]] = combdf
}
df_allalberta12s01 = do.call("rbind", listdf)

# Process COI

df_allalbertacoi = subset(df_allalberta02, grepl("COI", df_allalberta02$Method))
listdf = list()
i = 0 
for(method in unique(df_allalbertacoi$Method)){
  i = i + 1 
  sub = subset(df_allalbertacoi, Method == method)
  combdf = merge(list_truedf$coi, sub, by = "Species", all = TRUE)
  combdf = combdf %>%  fill(Method)
  listdf[[i]] = combdf
}
df_allalbertacoi01 = do.call("rbind", listdf)

# Combined both marker_gene

df_allalbertacomb01 = rbind(df_allalberta12s01, df_allalbertacoi01)

df_allalbertacomb01$Prediction_classification = df_allalbertacomb01$Prediction_classification %>% 
  replace_na("False positive")
df_allalbertacomb01$Prediction_classification = ifelse(is.na(df_allalbertacomb01$Total),
                                                       "False negative",
                                                       df_allalbertacomb01$Prediction_classification)
df_allalbertacomb02 = df_allalbertacomb01 %>% 
  group_by(Method, Prediction_classification) %>%
  summarise(Count = n())

df_allalbertacomb03 = dcast(df_allalbertacomb02,
                            Method ~ Prediction_classification, value.var = "Count", fill = 0)

df_allalbertacomb03$Precision = df_allalbertacomb03$`True positive` /
  (df_allalbertacomb03$`True positive` + df_allalbertacomb03$`False positive`)

df_allalbertacomb03$Recall = df_allalbertacomb03$`True positive` /
  (df_allalbertacomb03$`True positive` + df_allalbertacomb03$`False negative`)

df_allalbertacomb03$Sample = "All Alberta fish species"
df_allalbertacomb01$Sample = "All Alberta fish species"

############################ MIX ALBERTA ##################################

# Reusing the lists of Alberta fish included in the mock community generated previously (All Alberta) for both marker genes and adding the other fishes

otherfish = read.table("/mock_communities/names_otherfishmix.txt", sep = ";")
names(otherfish) = c("ID", "Species")
otherfish$Species = gsub("_", " ", otherfish$Species)
otherfish$Prediction_classification = "True positive"
otherfish = select(otherfish, c("Species", "Prediction_classification"))

true_other12s = rbind(otherfish, list_truedf$`12s`)
true_othercoi = rbind(otherfish, list_truedf$coi)


df_mix01 = subset(df01, Sample == "All Alberta and other species of fish" & 
                    !grepl("<", df01$Species))

df_mix01$Species = ifelse(df_mix01$Species == "Non-Alberta", 
                          df_mix01$original_specie, df_mix01$Species)

df_mix02 = df_mix01 %>% 
  group_by(Method, Species) %>%
  summarize(Total = n())

# Process 12S

df_mix12s = subset(df_mix02, grepl("12S", df_mix02$Method))

listdf = list()
i = 0 
for(method in unique(df_mix12s$Method)){
  i = i + 1 
  sub = subset(df_mix12s, Method == method)
  combdf = merge(true_other12s, sub, by = "Species", all = TRUE)
  combdf = combdf %>%  fill(Method)
  listdf[[i]] = combdf
}
df_mix12s01 = do.call("rbind", listdf)

# Process COI

df_mixcoi = subset(df_mix02, grepl("COI", df_mix02$Method))

listdf = list()
i = 0 
for(method in unique(df_mixcoi$Method)){
  i = i + 1 
  sub = subset(df_mixcoi, Method == method)
  combdf = merge(true_othercoi, sub, by = "Species", all = TRUE)
  combdf = combdf %>%  fill(Method)
  listdf[[i]] = combdf
}
df_mixcoi01 = do.call("rbind", listdf)


#### Verify is the specie is present in the different databases ####

# Load list of species in the different databases 

species_12scurated = list_truedf$`12s`
species_12sbarque = read.table("/database/header_barque_12S.txt", sep = "_")
species_12sbarque$Species = paste(species_12sbarque$V2, species_12sbarque$V3, sep = " ")
species_12sbarque$Species = gsub("Phoxinus neogaeus","Chrosomus neogaeus", species_12sbarque$Species)
species_12sbarque$Species = gsub("Phoxinus eos","Chrosomus eos", species_12sbarque$Species)
species_12sbarque[nrow(species_12sbarque) + 1, ] = c("Cottidae","Cottus", "bairdii","Cottus bairdii:cognatus")
species_12sbarque[nrow(species_12sbarque) + 1, ] = c("Salmonidae","Salvelinus", "alpinus","Salvelinus alpinus:malma")


species_12srdp = read.table("/database/header_rdp_12S.txt", sep = ";")
species_12srdp$V9 = gsub("_", " ", species_12srdp$V9)
species_12srdp[nrow(species_12srdp) + 1, ] = c("1","2", "3", "4", "5", "6", "7", "8","Cottus bairdii:cognatus")
species_12srdp[nrow(species_12srdp) + 1, ] = c("1","2", "3", "4", "5", "6", "7", "8","Salvelinus alpinus:malma")

species_coicurated = list_truedf$coi
species_coibarque = read.table("/database/header_bold_coi_barque.txt", sep = "_")
species_coibarque$Species = paste(species_coibarque$V2, species_coibarque$V3, sep = " ")
species_coibarque$Species = gsub("Hemichromis bimaculatus", "Rubricatochromis bimaculatus", species_coibarque$Species)
species_coibarque[nrow(species_coibarque) + 1, ] = c("1","2", "3","Coregonus artedi:zenithicus")

species_coirdp = read.table("/header_rdp_coi.txt", sep = ";")
species_coirdp$V9 = gsub("_", " ", species_coirdp$V9)
species_coirdp$V9 = gsub("Hemichromis bimaculatus", "Rubricatochromis bimaculatus", species_coirdp$V9)
species_coirdp[nrow(species_coirdp) + 1, ] = c("1","2", "3", "4", "5", "6", "7", "8", "Coregonus artedi:zenithicus")

# Cross reference with list of species in db 

df_12scurate = subset(df_mix12s01, grepl("Curated", df_mix12s01$Method))
df_12scurate$Prediction_classification = ifelse(!df_12scurate$Species %in% species_12scurated$Species, 
                                                "Not in database", df_12scurate$Prediction_classification)

df_coicurate = subset(df_mixcoi01, grepl("Curated", df_mixcoi01$Method))
df_coicurate$Prediction_classification = ifelse(!df_coicurate$Species %in% species_coicurated$Species, 
                                                "Not in database", df_coicurate$Prediction_classification)

df_12srdp = subset(df_mix12s01, Method == "12S - DADA2 - RDP - Non-curated")
df_12srdp$Prediction_classification = ifelse(!df_12srdp$Species %in% species_12srdp$V9, 
                                             "Not in database", df_12srdp$Prediction_classification)

df_12sbarque = subset(df_mix12s01, Method == "12S - Barque - VSEARCH - Non-curated")
df_12sbarque$Prediction_classification = ifelse(!df_12sbarque$Species %in% species_12sbarque$Species , 
                                                "Not in database", df_12sbarque$Prediction_classification)

df_coirdp = subset(df_mixcoi01, Method == "COI - DADA2 - RDP - Non-curated")
df_coirdp$Prediction_classification = ifelse(!df_coirdp$Species %in% species_coirdp$V9, 
                                             "Not in database", df_coirdp$Prediction_classification)

df_coibarque = subset(df_mixcoi01, Method == "COI - Barque - VSEARCH - Non-curated")
df_coibarque$Prediction_classification = ifelse(!df_coibarque$Species %in% species_coibarque$Species , 
                                                "Not in database", df_coibarque$Prediction_classification)

# Remove all species that are not in database. Any NA left are TRUE false negative 

df_mixall01 = rbind(df_12scurate, df_12srdp, df_12sbarque, df_coicurate, df_coirdp, df_coibarque)

df_mixall01$Prediction_classification = df_mixall01$Prediction_classification %>% 
  replace_na("False positive")

df_mixall02 = subset(df_mixall01, Prediction_classification != "Not in database")

df_mixall02$Prediction_classification = ifelse(is.na(df_mixall02$Total),
                                               "False negative",
                                               df_mixall02$Prediction_classification)
df_mixall03 = df_mixall02 %>% 
  group_by(Method, Prediction_classification) %>%
  summarise(Count = n())

df_mixall04 = dcast(df_mixall03,
                    Method ~ Prediction_classification, value.var = "Count", fill = 0)

df_mixall04$Precision = df_mixall04$`True positive` /
  (df_mixall04$`True positive` + df_mixall04$`False positive`)

df_mixall04$Recall = df_mixall04$`True positive` /
  (df_mixall04$`True positive` + df_mixall04$`False negative`)

df_mixall04$Sample = "All Alberta and other species of fish"
df_mixall02$Sample = "All Alberta and other species of fish"

##############################################################################
########################## Generated graphs ##################################  
##############################################################################

# Empirical data 

df_empmelt = melt(select(df_empirical, c("Method", "Precision","Sample")))

df_empmelt$Method = factor(df_empmelt$Method, levels = c("12S - Barque - VSEARCH - Curated", 
                                                         "12S - Barque - VSEARCH - Non-curated", 
                                                         "COI - Barque - VSEARCH - Curated", 
                                                         "COI - Barque - VSEARCH - Non-curated", 
                                                         "12S - DADA2 - RDP - Curated", 
                                                         "12S - DADA2 - RDP - Non-curated", 
                                                         "COI - DADA2 - RDP - Curated",         
                                                         "COI - DADA2 - RDP - Non-curated"))

# df_empmelt$Method = gsub("Barque - ", "Barque<br>", df_empmelt$Method)
# df_empmelt$Method = gsub("DADA2 - ", "DADA2<br>", df_empmelt$Method)

df_empmelt$Method = gsub(" - ", "<br>", df_empmelt$Method)

# In-Silico communities 

df_mock= rbind(df_mixall04, df_allalbertacomb03, df_low05, df_blindsplot05)

df_mockmelt = melt(select(df_mock, c("Method", "Recall", "Precision","Sample")))

df_mockmelt$Method = factor(df_mockmelt$Method, levels = c("12S - Barque - VSEARCH - Curated", 
                                                           "12S - Barque - VSEARCH - Non-curated", 
                                                           "COI - Barque - VSEARCH - Curated", 
                                                           "COI - Barque - VSEARCH - Non-curated", 
                                                           "12S - DADA2 - RDP - Curated", 
                                                           "12S - DADA2 - RDP - Non-curated", 
                                                           "COI - DADA2 - RDP - Curated",         
                                                           "COI - DADA2 - RDP - Non-curated"))

# df_mockmelt$Method = gsub("Barque - ", "Barque<br>", df_mockmelt$Method)
# df_mockmelt$Method = gsub("DADA2 - ", "DADA2<br>", df_mockmelt$Method)
df_mockmelt$Method = gsub(" - ", "<br>", df_mockmelt$Method)


dfmockrecall = subset(df_mockmelt, variable == "Recall")
dfmockprecision = subset(df_mockmelt, variable == "Precision")

plot_mock_recall = ggplot(dfmockrecall, aes(x = Method, y = value, color= Method, fill = Method)) + 
  geom_boxplot(alpha = 0.2, width = 0.3, linewidth = 0.3) + 
  geom_jitter(alpha = 0.8, size = 1) +
  custom_theme() + 
  theme(panel.border = element_blank(), 
        panel.spacing = unit(15, "pt"),
        axis.line.x = element_line(linewidth = 0.2),         
        axis.line.y = element_blank(), 
        legend.position = "none",
        axis.ticks.y = element_blank(),
        #axis.text.x = element_markdown(angle = 90, hjust = 0.95,  vjust = 0.5),
        axis.text.x = element_markdown(size = 6),      
        axis.title.x = element_blank()
        
  ) + 
  labs(y = "Recall") + 
  scale_color_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                                "#006400", "#70af2f", "orange",  "#592321"), name ="Method") + 
  
  scale_fill_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                               "#006400", "#70af2f", "orange",  "#592321"), name ="Method") 

plot_mock_precision = ggplot(dfmockprecision, aes(x = Method, y = value, color= Method, fill = Method)) + 
  geom_boxplot(alpha = 0.2, width = 0.3, linewidth = 0.3) + 
  geom_jitter(alpha = 0.8, size = 1) +
  custom_theme() + 
  theme(panel.border = element_blank(), 
        panel.spacing = unit(15, "pt"),
        axis.line.x = element_line(linewidth = 0.2),         
        axis.line.y = element_blank(), 
        legend.position = "none",
        axis.ticks.y = element_blank(),
        #axis.text.x = element_markdown(angle = 90, hjust = 0.95,  vjust = 0.5),
        axis.text.x = element_markdown(size = 6), 
        axis.title.x = element_blank()
        
  ) + 
  labs(y = "Precision") + 
  scale_color_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                                "#006400", "#70af2f", "orange",  "#592321"), name ="Method") + 
  
  scale_fill_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                               "#006400", "#70af2f", "orange",  "#592321"), name ="Method") 

plot_empirial = ggplot(df_empmelt, aes(x = Method, y = value, color = Method, fill = Method)) + 
  geom_boxplot(alpha = 0.2, width = 0.3, linewidth = 0.3) + 
  geom_jitter(alpha = 0.8, size = 1) +
  custom_theme() + 
  theme(panel.border = element_blank(), 
        axis.line.y = element_blank(), 
        axis.line.x = element_line(linewidth = 0.2),        
        legend.position = "none",
        axis.ticks.y = element_blank(),
        #axis.text.x = element_markdown(angle = 90, hjust = 0.95, vjust = 0.5), 
        axis.text.x = element_markdown(size = 6),        
        axis.title.x = element_blank(),
        
  ) + 
  labs(y = "Precision") + 
  scale_color_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                                "#006400", "#70af2f", "orange",  "#592321"), name ="Method") + 
  
  scale_fill_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                               "#006400", "#70af2f", "orange",  "#592321"), name ="Method") 

# Plot of proportion of assigned as True Alberta fish, Non-fish, Non-Alberta fish 

df_empproportion$Method = factor(df_empproportion$Method, levels = c("12S - Barque - VSEARCH - Curated", 
                                                                     "12S - Barque - VSEARCH - Non-curated", 
                                                                     "COI - Barque - VSEARCH - Curated", 
                                                                     "COI - Barque - VSEARCH - Non-curated", 
                                                                     "12S - DADA2 - RDP - Curated", 
                                                                     "12S - DADA2 - RDP - Non-curated", 
                                                                     "COI - DADA2 - RDP - Curated",         
                                                                     "COI - DADA2 - RDP - Non-curated"))


df_empproportion$Method = gsub("Barque - ", "Barque<br>", df_empproportion$Method)
df_empproportion$Method = gsub("DADA2 - ", "DADA2<br>", df_empproportion$Method)



plot_proportion = ggbarplot(df_empproportion, color = "Method", 
                            fill = "Method", y = "Proportion", 
                            x = "Group_taxonomy", alpha = 0.2, add = c("mean_se"),
                            position = position_dodge(width = 0.9, preserve = "single")) + 
  custom_theme() + 
  theme(panel.border = element_blank(), 
        axis.ticks = element_blank(),
        axis.line.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.2),   
        axis.text.y = element_text(hjust = 1),
        axis.text.x = element_markdown(),
        axis.title = element_blank(),
        legend.key.spacing.y = unit(0.1, "cm"),
        legend.text = element_markdown(size = 6),
        legend.position = "right", 
        legend.key.size = unit(0.5, 'cm')) +
  scale_y_continuous(expand = c(0,0), limits = c(0,105)) + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_color_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                                "#006400", "#70af2f", "orange",  "#592321"), name ="Method") + 
  
  scale_fill_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                               "#006400", "#70af2f", "orange",  "#592321"), name ="Method") 


plot = ggarrange(plot_mock_precision,plot_empirial,plot_mock_recall,plot_proportion, 
                 labels = c("A ", "B ", "C ", "D "), nrow = 2, ncol = 2, align = "h")

ggsave("/outputs/figures/all_recall_precision.pdf", plot, height = 5, width = 8.5, units = "in")  

all_df = rbind(dfmockrecall, dfmockprecision, df_empmelt)

write.csv(all_df, "/outputs/recall_precision.csv", quote = FALSE, row.names = FALSE)


###############################################################################
########################## Statistical tests ##################################  
###############################################################################

###### Precision empirical data set ###### 

df_empmelt$Method = gsub("<br>", " - ", df_empmelt$Method)

df_empmelt2 = df_empmelt %>%
  separate(col = Method, into = c("Marker_gene", "Workflow", "Classified", "Database"), sep = " - ")

# Compare curated VS non-curated 
ktest_database = kruskal.test(value ~ Database, data = df_empmelt2)

# Compared curated 12S vs curated COI
curated = subset(df_empmelt2, Database == "Curated")
kruskal.test(value ~ Marker_gene, data = curated)
curated_av_permarker_gene = curated %>%
  group_by(Marker_gene) %>%
  summarise(Average = mean(value), 
            sd = sd(value, na.rm = TRUE))

# Compared non-curated 12S vs curated COI
noncur = subset(df_empmelt2, Database == "Non-curated")
kruskal.test(value ~ Marker_gene, data = noncur)
noncur_av_permarker_gene = noncur %>%
  group_by(Marker_gene) %>%
  summarise(Average = mean(value), 
            sd = sd(value, na.rm = TRUE))


###### Recall ###### 

dfmockrecall$Method = gsub("<br>", " - ", dfmockrecall$Method)
df_recall = dfmockrecall %>%
  separate(col = Method, into = c("Marker_gene", "Workflow", "Classified", "Database"), sep = " - ")

ktest_database = kruskal.test(value ~ Database, data = df_recall)