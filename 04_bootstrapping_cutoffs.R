library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggtext)
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
      #anel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
      #panel.grid.minor = element_blank(), 
      #panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid', colour = "grey")
      #panel.grid.major = element_blank()
    )
}

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
                                    "Confidence score<br>below threhold", 
                                    "Confidence score<br>above threhold")

# Remove certain samples we are not using 
df01 = subset(df01, !grepl("richness|distribution", df01$Sample))

rate_confidence = df01 %>%
  group_by(Sample, Method, Confidence_categories) %>%
  summarise(Total_pergroup = n())

rate_confidence = rate_confidence %>% 
  group_by(Sample, Method) %>%
  mutate(Rate = Total_pergroup/sum(Total_pergroup)*100)

rate_confidence$Method = factor(rate_confidence$Method, levels = c("12S - Barque - VSEARCH - Curated", 
                                                                   "12S - Barque - VSEARCH - Non-curated", 
                                                                   "COI - Barque - VSEARCH - Curated", 
                                                                   "COI - Barque - VSEARCH - Non-curated", 
                                                                   "12S - DADA2 - RDP - Curated", 
                                                                   "12S - DADA2 - RDP - Non-curated", 
                                                                   "COI - DADA2 - RDP - Curated",         
                                                                   "COI - DADA2 - RDP - Non-curated"))

df_emp = subset(rate_confidence, grepl("00", rate_confidence$Sample))
df_mock = subset(rate_confidence, !grepl("00", rate_confidence$Sample))


plot_mock = ggbarplot(df_mock, color = "Method", 
                      fill = "Method", y = "Rate", 
                      x = "Confidence_categories", alpha = 0.2, add = c("mean_se"),
                      position = position_dodge(width = 0.8, preserve = "single")) + 
  custom_theme() + 
  ggtitle("BREAK<br>") + 
  theme(panel.border = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_markdown(hjust = 0, size = 12, color = "white"),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(),
        legend.position = "right", 
        legend.key.size = unit(0.5, 'cm')) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_color_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                                "#006400", "#70af2f", "orange",  "#592321"), name ="Method") + 
  
  scale_fill_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc",
                               "#006400", "#70af2f", "orange",  "#592321"), name ="Method") 


plot_empiricaldata = ggbarplot(df_emp, color = "Method", 
                               fill = "Method", y = "Rate", 
                               x = "Confidence_categories", alpha = 0.2, add = c("mean_se"),
                               position = position_dodge(width = 0.8, preserve = "single")) + 
  custom_theme() + 
  ggtitle("BREAK<br>") + 
  theme(panel.border = element_blank(), 
        plot.title = element_markdown(hjust = 0, size = 12, color = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_markdown(),
        legend.position = "right", 
        legend.key.size = unit(0.5, 'cm')) +
  scale_y_continuous(expand = c(0,1)) + 
  scale_x_discrete(expand = c(0,0.5)) + 
  scale_color_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                                "#006400", "#70af2f", "orange",  "#592321"), name ="Method") + 
  
  scale_fill_manual(values = c("#063971", "#48c9b0", "#d9b3ff", "#ba64fc", 
                               "#006400", "#70af2f", "orange",  "#592321"), name ="Method") 

plot = ggarrange(plot_mock, plot_empiricaldata, common.legend = TRUE, legend = "right", labels = c("A", "B"))

ggsave("/figures/allsamples_confidencelevels.pdf", plot, height = 2.5, width = 8, units = "in")  
write.csv(rate_confidence, "/outputs/allsamples_rateconfidence.csv", quote = FALSE, row.names = FALSE)