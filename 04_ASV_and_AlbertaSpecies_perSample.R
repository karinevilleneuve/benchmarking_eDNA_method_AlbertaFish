
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

# Load list of Alberta fresh water fishes for common name  
alberta_fish = read.csv("/database/fishbase_alberta_with_missing.csv", header = TRUE, check.names = FALSE)
# Create list to ID fish not from Alberta
list_alberta_fish = unique(alberta_fish$Species)
extra_species = c("Oncorhynchus mykiss x Salmo salar", "Cottus bairdii:cognatus", 
                  "Salvelinus alpinus:malma", "Coregonus artedi:zenithicus")
list_alberta_species = c(list_alberta_fish, extra_species)


############################################
############# Load csv files ###############  
############################################

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

df01$Species = gsub("Hemichromis bimaculatus", "Rubricatochromis bimaculatus", df01$Species)

df_02 = subset(df01, grepl("00", df01$Sample) &
                 !grepl("<", df01$Species))

df_02$Method = factor(df_02$Method, levels = c("12S - Barque - VSEARCH - Curated", 
                                               "12S - Barque - VSEARCH - Non-curated", 
                                               "COI - Barque - VSEARCH - Curated", 
                                               "COI - Barque - VSEARCH - Non-curated", 
                                               "12S - DADA2 - RDP - Curated", 
                                               "12S - DADA2 - RDP - Non-curated", 
                                               "COI - DADA2 - RDP - Curated",         
                                               "COI - DADA2 - RDP - Non-curated"))

df_02$Method = gsub(" - ", "<br>", df_02$Method)

# Count of unique species 

df_onlyalberta = subset(df_02, !grepl("-", df_02$Species) &
                          Abundance > 0)

df_specie_count = df_onlyalberta %>%
  group_by(Method, Sample) %>%
  summarise(Count_unique_species = n_distinct(Species))

unique_species = ggplot(df_specie_count, aes(x = Method, y = Count_unique_species, color = Method, fill = Method)) + 
  geom_boxplot(alpha = 0.2, width = 0.3, linewidth = 0.2) + 
  geom_jitter(alpha = 0.8, size = 1) +
  custom_theme() + 
  theme(panel.border = element_blank(), 
        axis.line.y = element_blank(), 
        axis.line.x = element_line(linewidth = 0.2),        
        legend.position = "none",
        axis.ticks.y = element_blank(),
        #axis.text.x = element_markdown(angle = 90, hjust = 0.95, vjust = 0.5), 
        axis.text.x = element_markdown(size = 6),        
        axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)) + 
  labs(y = "Unique species of Alberta fish") + 
  scale_color_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                                "#006400", "#70af2f", "orange",  "#592321"), name ="Method") + 
  
  scale_fill_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                               "#006400", "#70af2f", "orange",  "#592321"), name ="Method") 

############################################
########### TOTAL NUMBER ASV ###############  
############################################

df_asvtotal = df_02 %>%
  group_by(Method, Sample) %>%
  summarise(Total_asv = sum(Abundance))

plot_total = ggplot(df_asvtotal, aes(x = Method, y = Total_asv, color = Method, fill = Method)) + 
  geom_boxplot(alpha = 0.2, width = 0.3, linewidth = 0.2) + 
  geom_jitter(alpha = 0.8, size = 1) +
  custom_theme() + 
  theme(panel.border = element_blank(), 
        axis.line.y = element_blank(), 
        axis.line.x = element_line(linewidth = 0.2),        
        legend.position = "none",
        axis.ticks.y = element_blank(),
        #axis.text.x = element_markdown(angle = 90, hjust = 0.95, vjust = 0.5), 
        axis.text.x = element_markdown(size = 6),        
        axis.title.x = element_blank()
        
  ) + 
  labs(y = "Number of ASV") + 
  scale_color_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                                "#006400", "#70af2f", "orange",  "#592321"), name ="Method") + 
  
  scale_fill_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                               "#006400", "#70af2f", "orange",  "#592321"), name ="Method") 

plot = ggarrange(plot_total, unique_species, labels = c("A", "B"))

ggsave("/figures/asv_per_sample.pdf", plot, height = 2.5, width = 8.5, units = "in")


###############################################################################
########################## Statistical tests ##################################  
###############################################################################
library(rstatix)

df_asvtotal$Method = gsub("<br>", " - ", df_asvtotal$Method)

df_asvtotal2 = df_asvtotal %>%
  separate(col = Method, into = c("Marker_gene", "Workflow", "Classified", "Database"), sep = " - ")

# Comparing marker genes - Total number of ASV 
kruskal.test(Total_asv ~ Marker_gene, data = df_asvtotal2)

df_asvmean = df_asvtotal2 %>%
  group_by(Marker_gene) %>%
  summarise(Average = mean(Total_asv),
            std = sd(Total_asv, na.rm = TRUE))

m12s = subset(df_asvtotal2, Marker_gene == "12S") 
m12s$Method = paste(m12s$Workflow, m12s$Database)

kruskal.test(Total_asv ~ Method, data = m12s)

mcoi = subset(df_asvtotal2, Marker_gene == "COI") 
mcoi$Method = paste(mcoi$Workflow, mcoi$Database)

kruskal.test(Total_asv ~ Method, data = mcoi)

mcoi$method = paste(mcoi$Marker_gene, mcoi$Workflow, mcoi$Classified, mcoi$Database)

mcoi %>% 
  wilcox_test(Total_asv ~ method, p.adjust.method = "bonferroni") 

kruskal.test(Total_asv ~ Classified, data = mcoi)

mcoi %>% group_by(method) %>%
  summarise(average = mean(Total_asv),
            std = sd(Total_asv, na.rm = TRUE))

# PResence absence heatmap 

df_species = df_onlyalberta %>%
  group_by(Method, Species) %>%
  summarise(Count = n())

df_species$ok = 1 

df_species = df_species %>%
  mutate(Species = str_replace(Species, "(.*)","*\\1*")) # Adding asterisk at beginning and end of every taxa

plot_heatmap = ggplot(df_species, aes(y = Species, x = Method, fill = ok)) + 
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  custom_theme() + 
  theme(panel.border = element_blank(), 
        axis.line.y = element_blank(), 
        axis.line.x = element_line(linewidth = 0.2),        
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_markdown(hjust = 1),
        axis.text.x = element_markdown(size = 6),        
        axis.title.x = element_blank()) + 
  scale_y_discrete(limits=rev)

ggsave("/figures/species_per_method.pdf", plot_heatmap, height = 4, width = 5, units = "in")