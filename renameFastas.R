## Extract genes

# Packages
library(phytools)
library(stringr)

# Note to future Jun
# This file is usually run AFTER assembly but BEFORE annotations are extracted
# In this case, because of the suffix "bwa", the renaming did not succeed
# I am writing a stop-gap measure (right at the bottom) to fix it at the alignment stage
# In many ways, fixing the names later right than earlier makes more sense anyway

data_dir <- "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data"
#chloroplast_dir <- file.path(data_dir, "chloroplast_bwa_assembled")
#ribosome_dir <- file.path(data_dir, "ribosome_bwa_assembled")
#output_dir <- "~/Desktop/pep_130618"

chloroplastData <- read.dna(file.path(data_dir, "concatAlign_160618.phy"), "interleaved")
ncol(chloroplastData) # 63265

# Calculate number of parsimony informative sites
parsInfo <- vector()
for(i in 1:ncol(chloroplastData)){
  # for each aligned site
  count <- table(as.character(chloroplastData[,i]))
  count_subset <- count[! names(count) %in% "n"] # remove Ns
  count_parsInfo <- count_subset[which(count_subset > 2)]
  parsInfo[i] <- ifelse(length(count_parsInfo) > 1, 1, 0)
}
sum(parsInfo) / length(parsInfo) * 100
# 9.4% parsimony informative sites

# R =	A or G
# Y	= C or T
# S	= G or C
# W	= A or T
# K	= G or T
# M	= A or C
# B	= C or G or T
# D	= A or G or T
# H	= A or C or T
# V	= A or C or G
# N	any base or -	gap


# Accession data
accessionData <- read.csv("/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/data/PeperomiaAccessionData.csv")
accessionData$label <- paste(accessionData$SampleID, accessionData$Genus, accessionData$Species, accessionData$Geography,sep = "_")
accessionData$label <- gsub(accessionData$label, pattern = " ", replacement = "_")

# Clean up tip names
chloroplastLabels <- gsub(rownames(chloroplastData), pattern = "_bwa", replacement = "")
rownames(chloroplastData) <- accessionData$label[match(chloroplastLabels, accessionData$SampleID)]
rownames(chloroplastData)[is.na(rownames(chloroplastData))] <- "Piper_kadsura"

write.dna(chloroplastData, file.path(data_dir, "concatAlign_160618_clean.fasta"), "fasta")

# # Rename chloroplast files
# chloroplast_files <- list.files(chloroplast_dir)
# chloroplast_id <- str_split_fixed(chloroplast_files, pattern = "_", n = 2)[,1]
# 
# for(i in 1:length(chloroplast_files)){
#   temp <- read.dna(file.path(chloroplast_dir, chloroplast_files[i]), "fasta")
#   rownames(temp) <- accessionData$label[match(chloroplast_id[i], accessionData$SampleID)]
#   write.dna(temp, file.path(output_dir, chloroplast_files[i]), "fasta")
# }
# 
# ribosome_files <- list.files(ribosome_dir)
# ribosome_id <- str_split_fixed(ribosome_files, pattern = "_", n = 2)[,1]
# 
# for(i in 1:length(chloroplast_files)){
#   temp <- read.dna(file.path(ribosome_dir, ribosome_files[i]), "fasta")
#   rownames(temp) <- accessionData$label[match(ribosome_id[i], accessionData$SampleID)]
#   write.dna(temp, file.path(output_dir, ribosome_files[i]), "fasta")
# }
# 
# head(accessionData)
# peptaxa <- subset(accessionData, Genus == "Peperomia" & FinalLibrary.Y.N. == "Y")
# sort(unique(as.vector(peptaxa$Species)))
