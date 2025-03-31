library(dplyr)
library(stringr)
library(tidyr)
library(glue)
library(phyloseq)

###################################
####### Load fish list ############
###################################

# Load list of Alberta fresh water fishes for common name  
alberta_fish = read.csv("/database/fishbase_alberta_with_missing.csv", header = TRUE, check.names = FALSE)
# Create list to ID fish not from Alberta
list_alberta_fish = unique(alberta_fish$Species)

# Only keep column of interest for the common name 
alberta_fish_common_name = select(alberta_fish, c("Species", "Common name"))

# Load list of fish families downloaded from https://fishtreeoflife.org/taxonomy/
fish_fam = read.csv("/database/fish_families.csv", header = FALSE, check.names = FALSE)
list_fish_fam_raw = fish_fam$V1
# Add custom ID to list
extra_fam = c("Percentage ID < 97%", "Specie bootstrap < 70", "Leuciscidae",
              "Oncorhynchus_mykiss_x_Salmo_salar")
list_fish_fam = c(list_fish_fam_raw, extra_fam)

###################################
##### BARQUE Curated database #####
###################################

path_cur = "/barque_12s/results_12scurated"
files_annotation_cur = list.files(path = path_cur, pattern=".fasta.gz\\.") # Double backslash is important to esacape wildcard from period 
list_annotation_df_cur = list()
i = 0 

for (file in files_annotation_cur){
  i = i + 1
  sequencefile = gsub(".gz.12Scurated", "", file)
  sequence_raw = read.table(paste(path_cur, sequencefile, sep = "/"))
  seq_df = as.data.frame(matrix(sequence_raw$V1, ncol = 2, byrow = TRUE))
  seq_df$V1 = gsub(">", "", seq_df$V1)
  names(seq_df)  = c("otu_info", "sequence")
  annot_raw = read.table(file = paste(path_cur, file, sep = "/"))
  annotation = select(annot_raw, c("V1","V2", "V3"))
  names(annotation) = c("otu_info","classification", "Percentage ID")
  annotation = separate(annotation, classification, 
                        c("ID","CellularOrgamism", "Superkingdom", "Kingdom", "Phylum", 
                          "Class", "Order", "Family", "Genus","Species"), sep = ";")
  annotation_seq = merge(annotation, seq_df, by = "otu_info")
  annotation_clean = select(annotation_seq, c("otu_info", "Family", "Genus", "Species", "Percentage ID", "sequence"))
  annotation_clean$Database = "Curated"
  split_name = str_split_1(file, "_")
  sample_name = split_name[1]
  annotation_clean$Sample = sample_name
  list_annotation_df_cur[[i]] = annotation_clean
}  
df_all_annotation_cur = do.call("rbind", list_annotation_df_cur)

###################################
### BARQUE Non-Curated database ###
###################################

# The non-curate database does not include the full taxonomic rank, so we will be merging the results with a custom generate table of ranks. 

taxrank = read.csv("/database/reference_taxonomic_ranks.csv", header=TRUE)

path_raw = "/barque_12s/results_12sraw"
files_annotation_raw = list.files(path = path_raw, pattern=".fasta.gz\\.") # Double backslash is important to esacape wildcard from period 
list_annotation_df_raw = list()
i = 0 

for (file in files_annotation_raw){
  i = i + 1
  sequencefile = gsub(".gz.12S", "", file)
  sequence_raw = read.table(paste(path_raw, sequencefile, sep = "/"))
  seq_df = as.data.frame(matrix(sequence_raw$V1, ncol = 2, byrow = TRUE))
  seq_df$V1 = gsub(">", "", seq_df$V1)
  names(seq_df)  = c("otu_info", "sequence")
  annot_raw = read.table(file = paste(path_raw, file, sep = "/"))
  annotation = select(annot_raw, c("V1","V2", "V3"))
  names(annotation) = c("otu_info","classification", "Percentage ID")
  annotation = separate(annotation, classification, 
                        c("BOLD_rank", "Genus", "Epithet"), sep = "_")
  annotation$Species = paste(annotation$Genus, annotation$Epithet, sep = "_")
  annotation$Species = gsub("Phoxinus_", "Chrosomus_", annotation$Species)
  annotation_rank = merge(annotation, taxrank, by = "Species", all.y = FALSE, all.x = TRUE, suffixes = c(".x", ""))
  annotation_seq = merge(annotation_rank, seq_df, by = "otu_info")
  annotation_clean = select(annotation_seq, c("otu_info", "Family", "Genus", "Species", "Percentage ID", "sequence"))
  annotation_clean$Database = "Non-curated"
  split_name = str_split_1(file, "_")
  sample_name = split_name[1]
  annotation_clean$Sample= sample_name
  list_annotation_df_raw[[i]] = annotation_clean
}  
df_all_annotation_raw = do.call("rbind", list_annotation_df_raw)


#######################################
#####  Combine barque data frames #####
#######################################

df01_allannotation = rbind(df_all_annotation_cur, df_all_annotation_raw)

# For each ASV keep only the top hits 
df02_top_hit = df01_allannotation %>% 
  group_by(Database, Sample, otu_info) %>% 
  arrange(desc(`Percentage ID`)) %>%
  filter(`Percentage ID` == max(`Percentage ID`))


df03_top_hit_nodup = df02_top_hit[!duplicated(df02_top_hit),]


# Remove single and double-ton 
df03_top_hit_nodup = separate(df03_top_hit_nodup, otu_info, 
                              c("trash1", "otuID", "trash2", "Abundance", "trash3"), sep = "_")
df03_top_hit_nodup$Abundance = as.integer(df03_top_hit_nodup$Abundance)
df04 = subset(df03_top_hit_nodup, Abundance > 2)

# If ID is below 97% change Species and Family to Percentage ID < 97%
df04$Species = ifelse(df04$`Percentage ID` < 97, "Percentage ID < 97%", df04$Species)
df04$Family = ifelse(df04$Species == "Percentage ID < 97%", "Percentage ID < 97%", df04$Family)

# For hybrid specie, replace family with the Species name
df04$Family = ifelse(grepl("_x_", df04$Species), df04$Species, df04$Family)
df04$Species = gsub ("_", " ", df04$Species)

#######################################
######### Assess multiple hits ########
#######################################

duplicate = df04 %>% 
  group_by(Database, Sample, otuID) %>%
  filter(n()>1)


# Delete all multiple hits
df05_nodup = df04 %>% 
  group_by(Database, Sample, otuID) %>%
  filter(n()<2)

# Calculate relative abundances 
df_all_barque = df05_nodup %>% 
  group_by(Database, Sample) %>% 
  mutate(Relative_abundance = Abundance / sum(Abundance)) %>%
  ungroup()

df_all_barque$Method = paste("Barque - VSEARCH - ", df_all_barque$Database, sep = "")

#######################################
############# RDP output ##############
#######################################

path = "/outputs"
files_dada2_df = as.data.frame(t(list.files(path, full.names = TRUE, pattern = "12s_dada2")))

database_option = c("cur", "raw")

list_dada2 = list()
i = 0 
# Parsing the phyloseq object generated throught DADA2 
for (ps_obj in grep("phylo", files_dada2_df, value = TRUE)){
  i = i + 1
  ps = readRDS(file = grep(database_option[i], ps_obj, value = TRUE))
  meltps = psmelt(ps)
  df_clean = select(meltps, c("Sample", "Abundance","Family", "Species"))
  df_clean = df_clean %>% 
    group_by(Sample) %>% 
    mutate(Relative_abundance = Abundance / sum(Abundance)) %>%
    ungroup()
  df_clean$Database = database_option[i]
  list_dada2[[i]] = df_clean
}

df_all_rdp = do.call("rbind", list_dada2)
df_all_rdp$Database = ifelse(df_all_rdp$Database == "raw", "Non-curated", "Curated") 
df_all_rdp$Method = paste("DADA2 - RDP - ", df_all_rdp$Database, sep = "") 

#######################################
####### Combine RDP and barque ########
#######################################

# Keep only common columns 
commoncol = intersect(names(df_all_rdp), names(df_all_barque))
df_all_barque_clean = df_all_barque %>% select(commoncol)
# Combine data frames 
df01_combined = rbind(df_all_rdp, df_all_barque_clean)
df01_combined$original_specie = df01_combined$Species

# If species is not in list of fish Family from FishBase replace specie name with Non-fish
df01_combined$Species = ifelse(df01_combined$Family %in% list_fish_fam, df01_combined$Species, "Non-fish")

# If species is not in list of Alberta fish from FishBase replace specie name with Non-Alberta
list_otherspecies = c("Non-fish", "Percentage ID < 97%", "Specie bootstrap < 70", "Salvelinus alpinus:malma", "Cottus bairdii:cognatus")
list_alberta_fish2 = c(list_alberta_fish, list_otherspecies)

df01_combined$Species = ifelse(df01_combined$Species %in% list_alberta_fish2, df01_combined$Species, "Non-Alberta") 

# Combine with data frame of common names 
df02_comnames = merge(df01_combined, alberta_fish_common_name, by = "Species", all.y = FALSE, all.x = TRUE)
df02_comnames$Sample = gsub("-mifish", "", df02_comnames$Sample)

# Save dataframe 
write.csv(df02_comnames,"/outputs/12s_allmethods.csv", quote = FALSE, row.names = FALSE)