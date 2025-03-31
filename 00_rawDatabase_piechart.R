
# Script used to generate stack bar charts representing composition of databases ----- ####

#### ----------------- Extract sequence header from fasta ------------------------------------- ####

# Header were extracted using shell command : grep ">" file.fasta | sed -e "s/>//" > headers.txt

#### ------------------ Load libraries  ------------------------------------------------------- ####

library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)
library(randomcoloR) 
library(ggpubr) 

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

#### ----------------- Load taxonomic rank reference  ----------------------------------------- ####

ref = read.csv("/database/reference_taxonomic_ranks.csv", header = TRUE)

#### ----------------- Load databases sequence header ----------------------------------------- ####

## RDP classifier
raw_rdp12s = read.table("/database/header_rdp_12S.txt", sep = ";")
raw_rdpcoi = read.table("/database/header_rdp_coi.txt", sep = ";")

## barque
raw_barque12s = read.table("/database/header_barque_12S.txt", sep = "_")
raw_barquecoi = read.table("/database/header_bold_coi_barque.txt",  sep = "_")

#### :::: RDP databases:::: ####

# Which Orders are found in Class Undef_Chordata
undef_chord = subset(raw_rdpcoi, V5 == "undef_Chordata")
sort(unique(undef_chord$V6))

# [1] "Ceratodontiformes" "Crocodylia"        "Testudines"   
# Ceratodontiformes belongs in Class Dipnoi
# Crocodylia and Testudines belong in Class Reptilia 

raw_rdpcoi$V5 = ifelse(raw_rdpcoi$V5 == "undef_Chordata" & raw_rdpcoi$V6 == "Ceratodontiformes", "Dipnoi", raw_rdpcoi$V5)

raw_rdpcoi$V5 = ifelse(raw_rdpcoi$V5 == "undef_Chordata", "Reptilia", raw_rdpcoi$V5)
raw_rdpcoi$V5 = gsub("Lepidosauria", "Reptilia", raw_rdpcoi$V5)   

list_rdp_db = list("COI RDP" = raw_rdpcoi, "12S RDP" = raw_rdp12s)
list_rdp_taxa = list()

i = 0
for (databases in list_rdp_db){
  i = i + 1
  names(databases) = c("ID", "Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  databases$Class = ifelse(databases$Phylum != "Chordata", "Non-Chordata", databases$Class)
  df_grouped = databases %>% 
    group_by(Class) %>% 
    summarise(n = n()) %>% 
    as.data.frame()
  df_grouped$Abundance = df_grouped$n / sum(df_grouped$n)
  df_grouped$database = names(list_rdp_db[i])
  #  df_ref = merge(df_grouped, ref, by = "Genus", all.y = FALSE, all.x = TRUE)
  list_rdp_taxa[[i]] = df_grouped
  names(list_rdp_taxa)[i] = names(list_rdp_db[i])
}

df_rdp_final = do.call("rbind", list_rdp_taxa)

#### :::: barque databases:::: ####

list_barque_db = list("12S barque" = raw_barque12s, "COI barque" = raw_barquecoi)
list_barque_taxa = list()

i = 0
for (databases in list_barque_db){
  i = i + 1
  names(databases) = c("Rank1", "Genus", "sort_species")
  databases$Species = paste(databases$Genus,"_",databases$sort_species, sep ="")
  databases$database = names(list_barque_db[i])
  list_barque_taxa[[i]] = databases
}


# Combine dataframe of both databases together 
df_all_barque = do.call("rbind", list_barque_taxa)

# Merge with reference database to get taxonomic ranks 
df_barque_mergeref = merge(df_all_barque, ref, by = "Species", all.x = TRUE, all.y = FALSE, suffixes = c(x = "_barque", y = ""))

# Get rid of column sort_species Kingdom and Superkingdom
df_barque_mergeref = select(df_barque_mergeref, -c("sort_species", "Superkingdom", "Kingdom"))

# Correcting simple misclassification and typos 
df_barque_mergeref$Genus_barque = gsub("Ambistoma", "Ambystoma", df_barque_mergeref$Genus_barque)

df_barque_mergeref$Species = gsub("Dendroica_coronata", "Setophaga_coronata", df_barque_mergeref$Species)
df_barque_mergeref$Genus_barque = gsub("Dendroica", "Setophaga", df_barque_mergeref$Genus_barque)

df_barque_mergeref$Species = gsub("Peronophythora_litchii","Phytophthora_litchii", df_barque_mergeref$Species)
df_barque_mergeref$Genus_barque = gsub("Peronophythora", "Phytophthora",df_barque_mergeref$Genus_barque)

df_barque_mergeref$Species = gsub("Catla_catla","Labeo_catla", df_barque_mergeref$Species)
df_barque_mergeref$Genus_barque = gsub("Catla", "Labeo",df_barque_mergeref$Genus_barque)

df_barque_mergeref$Genus_barque = gsub("Chelonodon","Chelonodontops",df_barque_mergeref$Genus_barque)
df_barque_mergeref$Genus_barque = gsub("Chelonodontopstops","Chelonodontops",df_barque_mergeref$Genus_barque)

df_barque_mergeref$Species = gsub("Paraneetroplus_synspilus","Vieja_melanurus", df_barque_mergeref$Species)
df_barque_mergeref$Genus_barque = gsub("Paraneetroplus", "Vieja",df_barque_mergeref$Genus_barque)

df_barque_mergeref$Genus_barque = gsub("Emydoidea", "Emys",df_barque_mergeref$Genus_barque)


# Use genus with known higher ranks to fill missing ranks 
dfbarque_fill =   df_barque_mergeref %>%
  group_by(Genus_barque) %>% 
  fill(Genus, .direction = 'downup') %>%
  fill(Family, .direction = 'downup') %>%
  fill(Order, .direction = 'downup') %>%
  fill(Class, .direction = 'downup') %>%
  fill(Phylum, .direction = 'downup') %>%
  ungroup()

dfbarque_fill$Rank1 = str_to_title(dfbarque_fill$Rank1)

# Looking at missing Phylum 

missing_phylym = dfbarque_fill %>%
  group_by(Genus_barque) %>% 
  filter(is.na(Phylum))

rankNAphylum = unique(missing_phylym$Rank1)
uniquephylum = unique(dfbarque_fill$Phylum)
rankNOTphylum = setdiff(rankNAphylum, uniquephylum)
sort(rankNOTphylum)

sub = subset(dfbarque_fill, Rank1 %in% rankNOTphylum) 
no_ranks = unique(sub[!complete.cases(sub), ]) # 29 entries 

# Copied the dataframe no_ranks into excel and added the missing rank information
# Load missing infos 

ranks = read.table("/database/update_missing_rank.txt", sep = "\t", header = TRUE)

filter_comm = match(dfbarque_fill$Species, ranks$Species)
# Fill in the missing ranks using information in the ranks data frame 

dfbarque_fill$Genus = ifelse(is.na(dfbarque_fill$Genus), ranks$Genus[filter_comm], dfbarque_fill$Genus)

dfbarque_fill$Family = ifelse(is.na(dfbarque_fill$Family), ranks$Family[filter_comm], dfbarque_fill$Family)

dfbarque_fill$Order = ifelse(is.na(dfbarque_fill$Order), ranks$Order[filter_comm], dfbarque_fill$Order)

dfbarque_fill$Class = ifelse(is.na(dfbarque_fill$Class), ranks$Class[filter_comm], dfbarque_fill$Class)

dfbarque_fill$Phylum = ifelse(is.na(dfbarque_fill$Phylum), ranks$Phylum[filter_comm], dfbarque_fill$Phylum)

dfbarque_fill$Phylum = ifelse(is.na(dfbarque_fill$Phylum), dfbarque_fill$Rank1, dfbarque_fill$Phylum)

missing_chords = dfbarque_fill %>%
  group_by(Genus_barque) %>% 
  filter(is.na(Class)) %>%
  filter(Phylum == "Chordata")

sort(unique(missing_chords$Genus_barque)) # More than 1000 ... (zut)

# For these taxon we will assign to them the class "Unavailable Class of chordates" 

dfbarque_fill$Class = ifelse(dfbarque_fill$Phylum != "Chordata", "Non-Chordata", dfbarque_fill$Class)
dfbarque_fill$Class = ifelse(is.na(dfbarque_fill$Class), "Unavailable Class of chordate", dfbarque_fill$Class)


# What are undef_Chordata
undef = subset(dfbarque_fill, Class == "undef_Chordata")
unique(undef$Order)
# [1] "Testudines" "Crocodylia"
# Which are actually part of Class Reptilia
dfbarque_fill$Class = gsub("undef_Chordata", "Reptilia", dfbarque_fill$Class)
# Lepidosauria are also actually of subclass of Reptilia 
dfbarque_fill$Class = gsub("Lepidosauria", "Reptilia", dfbarque_fill$Class)

# Prepare the dataframe for plotting 
df_barque_final = dfbarque_fill %>% 
  group_by(database, Class) %>% 
  summarise(n = n()) %>%
  mutate(Abundance = n /sum(n)) %>% 
  as.data.frame()

df_final = rbind(df_barque_final, df_rdp_final)
sort(unique(df_final$Class))

df_final$Class = gsub("Appendicularia_tunicates_class", "Appendicularia", df_final$Class)

n = length(unique(df_final$Class)) 

# generate a set of X unique colors corresponding to the number of unique taxa
palette = distinctColorPalette(n)

# Add species name to color palette
names(palette) = unique(df_final$Class)

# assign gray to specific categories. The same nomenclature can be use to manually change certain colors.

# Non-chordate and missing information 
palette[["Non-Chordata"]] = "#E1E1E1"
palette[["Unavailable Class of chordate"]] = "black"
# Classes if "Fish" 
palette[["Actinopteri"]] = "#85AEDF"
palette[["Coelacanthimorpha"]] =   "#111E6C"
palette[["Dipnoi"]] =  "#003152"
palette[["Hyperoartia"]] = "#7285A5"
palette[["Chondrichthyes"]] = "#4C516D"
palette[["Cladistia"]] = "#a8dfe8"
palette[["Myxini"]] =  "#4F97A3"
# Others 
palette[["Amphibia"]] = "#C1E1C1"
palette[["Appendicularia"]] = "#ffc2c5"
palette[["Ascidiacea"]] = "#DDA7DA"
palette[["Aves"]] = "#D6AB9C"
palette[["Reptilia"]] =  "#74D176"
palette[["Leptocardii"]] =  "#D98289"
palette[["Mammalia"]] =  "#DB9D53"
palette[["Thaliacea"]] =  "#A841E7"


# Ordering the legend in alphabetical order
legend_raw = unique(df_final$Class) #Extract legend as text
ordered_legend = sort(legend_raw) # order alphabetically
reordered_legend = fct_relevel(ordered_legend, "Non-Chordata", "Unavailable Class of chordate") # move Others to the beginning
final_legend = levels(reordered_legend) # Extract the levels in a new object

my_scale = scale_fill_manual(name = "Class", breaks = paste(final_legend), values = palette, na.translate=FALSE, drop=TRUE, limits = force)

plot = ggplot(df_final, aes(x = database, weight = Abundance, fill = fct_reorder(Class, Abundance,.desc=FALSE))) +
  geom_bar() + 
  custom_theme() +
  labs(y ='Relative abundances', x="Database", title = "") +
  theme(text = element_text(size = text_size),
        panel.spacing = unit(2, "lines"),
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(size=text_size, vjust=0.5),
        axis.text.x = element_text(size=text_size, hjust=0.6), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust =0.5, size=text_size),
        axis.title=element_text(size=text_size),
        legend.position = "none",
        legend.title = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        legend.key.size = unit(0.5, 'cm'),
        legend.margin = margin(), # pre-emptively set zero margins 
        strip.background = element_blank()) + 
  my_scale + # Load our color palette 
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, by = 1), expand = c(0,0), labels = c("0", "1")) + # Remove the white space 
  # Adjusting the legend, notably the number of rows and position
  guides(fill = guide_legend(title = "Class", title.position = "top", reverse=FALSE)) 

# Save figure
ggsave("/figures/databases_composition.pdf", plot = plot, height = 2, width = 8, units = "in")

#trying donut chart 
donut = ggplot(df_final, aes(x = "", y = Abundance, fill = fct_reorder(Class, Abundance,.desc=FALSE))) +
  geom_bar(stat="identity", color = "white", width = 1) +
  coord_polar(theta = "y") +
  #  custom_theme() +
  facet_grid(~ database) + 
  theme(panel.spacing = unit(0.5, "lines"),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 7),
        legend.position = "none",
        strip.background = element_blank()) + 
  my_scale 

ggsave("/figures/databases_composition_piechart.pdf", plot = donut, height = 2, width = 8, units = "in")

# Legend by group
palette[["Actinopteri"]] = "#85AEDF"
palette[["Coelacanthimorpha"]] =   "#111E6C"
palette[["Dipnoi"]] =  "#003152"
palette[["Hyperoartia"]] = "#7285A5"
palette[["Chondrichthyes"]] = "#4C516D"
palette[["Cladistia"]] = "#a8dfe8"
palette[["Myxini"]] =  "#4F97A3"

class_fish = c("Actinopteri", "Coelacanthimorpha", "Dipnoi", "Hyperoartia", "Chondrichthyes", "Cladistia", "Myxini")

df_fish = subset(df_final, Class %in% class_fish)

plot_fish = ggplot(df_fish, aes(x = database, weight = Abundance, fill = fct_reorder(Class, Abundance,.desc=FALSE))) +
  geom_bar() + 
  custom_theme() +
  theme(legend.position = "right",
        legend.title = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        legend.key.size = unit(0.5, 'cm'),
        legend.margin = margin(), # pre-emptively set zero margins 
        strip.background = element_blank()) + 
  my_scale + # Load our color palette 
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, by = 0.5), expand = c(0,0), labels = c("0", "0.5", "1")) + # Remove the white space 
  # Adjusting the legend, notably the number of rows and position
  guides(fill = guide_legend(title = "Classes - Fish", title.position = "top", reverse=FALSE))

leg = get_legend(plot_fish)
leg_plot = as_ggplot(leg)

# Save figure
ggsave("/figures/databases_composition_fishlegend.pdf", plot = leg_plot, height = 5, width = 7, units = "in")


df_others = subset(df_final, !Class %in% class_fish)

plot_others = ggplot(df_others, aes(x = database, weight = Abundance, fill = fct_reorder(Class, Abundance,.desc=FALSE))) +
  geom_bar() + 
  custom_theme() +
  theme(legend.position = "right",
        legend.title = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        legend.key.size = unit(0.5, 'cm'),
        legend.margin = margin(), # pre-emptively set zero margins 
        strip.background = element_blank()) + 
  my_scale + # Load our color palette 
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, by = 0.5), expand = c(0,0), labels = c("0", "0.5", "1")) + # Remove the white space 
  # Adjusting the legend, notably the number of rows and position
  guides(fill = guide_legend(title = "Classes - Non-fish", title.position = "top", reverse=FALSE))

leg = get_legend(plot_others)
leg_plot = as_ggplot(leg)

ggsave("/databases_composition_legendothers.pdf", plot = leg_plot, height = 5, width = 7, units = "in")

#Save table 

write.csv(df_final, "/outputs/database_composition.csv", quote = FALSE, row.names = FALSE)