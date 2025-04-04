# Script used to parse the annotation output from RDP and combine these results 
# with the phyloseq object generated by the DADA2 worklfow. 

library(phyloseq)

columns_rdp_out_cur = c("ASV", "rank0","Root", "rank00", "rank_bootstrap", "CellularOrganisms", "rank1", "cellularOrganisms_bootstrap", "Superkingdom", "rank2", "superkingdom_bootstrap", "Kingdom", "rank3", "kingdom_bootstrap", "Phylum", "rank4", "phylum_bootstrap", "Class", "rank5","class_bootstrap", "Order", "rank6", "order_bootstrap", "Family", "rank7", "family_bootstrap", "Genus", "rank8", "genus_bootstrap", "Species", "rank9", "species_bootstrap")

columns_rdp_out_raw = c("ASV", "rank_bootstrap", "CellularOrganisms", "rank1", "cellularOrganisms_bootstrap", "Superkingdom", "rank2", "superkingdom_bootstrap", "Kingdom", "rank3", "kingdom_bootstrap", "Phylum", "rank4", "phylum_bootstrap", "Class", "rank5","class_bootstrap", "Order", "rank6", "order_bootstrap", "Family", "rank7", "family_bootstrap", "Genus", "rank8", "genus_bootstrap", "Species", "rank9", "species_bootstrap")

path_raw = "/raw_output"
path_save = "/outputs"

# Create dataframe with the different files
dada2_out = as.data.frame(t(list.files(path_raw, pattern = "_annotation.out", full.names = TRUE)))
list_ps = list()
i = 0 

# :::::: COI marker gene :::::: # 
for (files in dada2_out){
  i = i + 1
  # Load RDP annotation output 
  annot_rdp_raw = read.table(files, sep = "\t")
  # Rename columns 
  if (ncol(annot_rdp_raw) < 30){
    names(annot_rdp_raw) = columns_rdp_out_raw} else {
      names(annot_rdp_raw) = columns_rdp_out_cur}
  # If Bootstrap specie is under 0.70 replace specie and family with Specie bootstrap < 70
  annot_rdp_raw$Species = ifelse(annot_rdp_raw$species_bootstrap < 0.70, 
                                 "Specie bootstrap < 70", annot_rdp_raw$Species)
  annot_rdp_raw$Family = ifelse(annot_rdp_raw$species_bootstrap < 0.70, 
                                "Specie bootstrap < 70", annot_rdp_raw$Family)  
  # Replace space with underscore
  annot_rdp_raw$Species = gsub("_", " ", annot_rdp_raw$Species) 
  # Removing unnecessary columns 
  annot_matrix = annot_rdp_raw[,!grepl("rank|bootstrap|Root", names(annot_rdp_raw))]
  # ASV as row names 
  row.names(annot_matrix) = annot_matrix$ASV
  # Use file name to load phyloseq .rds object 
  filename_raw = strsplit(basename(files), "/")
  rdsfilename = paste(path_raw, "/", gsub("_rdp.*", "", filename_raw), "_phylo.rds", sep = "")
  # Load RDS object (phyloseq) 
  rds_obj = readRDS(rdsfilename)
  # Merge phyloseq object with taxonomy table (RDP output)
  ps = merge_phyloseq(rds_obj, tax_table(as.matrix(annot_matrix)))
  # Name to save file 
  savefilename = paste(path_save, "/", gsub("_annotation.out", "_phylo.rds", filename_raw), sep ="")
  saveRDS(ps, savefilename)
}