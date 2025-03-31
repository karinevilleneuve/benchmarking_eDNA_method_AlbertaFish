library(dada2)
library(ShortRead)
library(Biostrings)
library(reticulate)
library(dplyr)
library(phyloseq)

path = "/raw_FASTQ/12s"

myfiles = list.files(path = path, pattern="fastq.gz")

fnFs = sort(list.files(path,pattern = "_R1_001.fastq.gz",full.names=TRUE))
fnRs = sort(list.files(path,pattern = "_R2_001.fastq.gz",full.names=TRUE))
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

################################################################################
#### ------------------------ Removing primers ---------------------------- ####
################################################################################

FWD = "CCGGTAAAACTCGTGCCAGC"
REV = "CATAGTGGGGTATCTAATCCCAGTTTG"

allOrients = function(primer){
  require(Biostrings)
  dna = DNAString(primer)
  orients = c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients,toString))
}

FWD.orients = allOrients(FWD)
REV.orients = allOrients(REV)

primerHits = function(primer,fn){
  nhits = vcountPattern(primer,sread(readFastq(fn)),fixed=FALSE)
  return(sum(nhits>0))
}

rbind(FWD.ForwardReads=sapply(FWD.orients,primerHits,fnFs[[1]]),
      FWD.ReverseReads=sapply(FWD.orients,primerHits,fnRs[[1]]),
      REV.ForwardReads=sapply(REV.orients,primerHits,fnFs[[1]]),
      REV.ReverseReads=sapply(REV.orients,primerHits,fnRs[[1]]))


cutadapt = "/usr/local/bin/cutadapt"
system2(cutadapt, args = '--version')

path.cut = file.path(path, "cutadapt") # create folder cutadapt in directory path_raw
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut = file.path(path.cut, basename(fnFs))
fnRs.cut = file.path(path.cut, basename(fnRs))

FWD.RC = dada2:::rc(FWD)
REV.RC = dada2:::rc(REV)

R1.flags = paste("-g", FWD, "-a", REV.RC)
R2.flags = paste("-G", REV, "-A", FWD.RC)

#Trim primers with cutadapt
for(i in seq_along(fnFs)){
  system2(cutadapt, arg = c(R1.flags, R2.flags, "-n", 2, "-m", 10, "-o", 
                            fnFs.cut[i],"-p",fnRs.cut[i],fnFs[i],fnRs[i]))
}

# Sanity check to confirm primers were removed 
rbind(FWD.ForwardReads=sapply(FWD.orients,primerHits,fnFs.cut[[1]]),
      FWD.ReverseReads=sapply(FWD.orients,primerHits,fnRs.cut[[1]]),
      REV.ForwardReads=sapply(REV.orients,primerHits,fnFs.cut[[1]]),
      REV.ReverseReads=sapply(REV.orients,primerHits,fnRs.cut[[1]]))

################################################################################
#### ------------------------- Filter and trim ---------------------------- ####
################################################################################

cutFs = sort(list.files(path.cut, pattern ="_R1_001.fastq.gz", full.names=TRUE))
cutRs = sort(list.files(path.cut, pattern ="_R2_001.fastq.gz", full.names=TRUE))


filtFs = file.path(path,"filtered", basename(cutFs))
filtRs = file.path(path,"filtered", basename(cutRs))

out = filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2,2),truncQ = 2,
                    minLen = 50, rm.phix = TRUE, compress = TRUE)

head(out)
tail(out)
exists = file.exists(filtFs) # All files have reads

errF = learnErrors(filtFs[exists], multithread = TRUE)
errR = learnErrors(filtRs[exists], multithread = TRUE)


derepFs = derepFastq(filtFs[exists],verbose=TRUE)
derepRs = derepFastq(filtRs[exists],verbose=TRUE)

dadaFs = dada(derepFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs = dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
names(mergers) = sample.names[exists]
seqtab = makeSequenceTable(mergers)
table(nchar(getSequences(seqtab))) # Most sequences are 253 bp long
#seqtab2 = seqtab[,nchar(colnames(seqtab)) %in% 200:300]

seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE,verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab) # 98.8% of reads remain

getN = function(x) sum(getUniques(x))
out_ex = row.names(out)[exists]

track = cbind(out,sapply(dadaFs,getN),sapply(dadaRs,getN),sapply(mergers,getN),
              rowSums(seqtab.nochim))

colnames(track) = c("input","filtered","denoisedF","denoisedR","merged","nochim")

ps_base = phyloseq(otu_table(seqtab.nochim,taxa_are_rows = FALSE))
dna = Biostrings::DNAStringSet(taxa_names(ps_base))
names(dna) = taxa_names(ps_base)
ps_base = merge_phyloseq(ps_base,dna)
ps_sub = prune_taxa(taxa_sums(ps_base) > 2, ps_base) # 2982 taxa remain

taxa_names(ps_sub)<-paste0("ASV",seq(ntaxa(ps_sub)))
ps_sub %>%
  refseq() %>%
  Biostrings::writeXStringSet("/raw_output/12s_dada2_biostring4annotation.fna", append = FALSE, compress = FALSE, compression_level = NA, format = "fasta")

saveRDS(ps_sub, file = "/raw_output/12s_dada2_phylo.rds")

################################################################################
#### ----------------------------- Classify ------------------------------- ####
################################################################################

# RDP classifer to assign taxonomy to 12S ASVs from 12S fish classifier https://github.com/terrimporter/12SfishClassifier/tree/v1.0.1
# java -Xmx1g -jar ~/rdp_classifier_2.14/dist/classifier.jar classify -t /12Sfishclassifier/mydata_trained/rRNAClassifier.properties -o /raw_output/12s_dada2_rdpraw_annotation.out -q /raw_output/12s_dada2_biostring4annotation.fna
#RDP classifier with curated database 
# java -Xmx1g -jar ~/rdp_classifier_2.14/dist/classifier.jar classify -t /curateddatabase/curated12s_db/12s_training_files/rRNAClassifier.properties -o /raw_output/12s_dada2_rdpcur_annotation.out -q /raw_output/12s_dada2_biostring4annotation.fna


