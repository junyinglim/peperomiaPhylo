## Extract genes

# Packages
library(phytools)
library(stringr)

data_dir <- "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data"
chloroplast_dir <- file.path(data_dir, "chloroplast_bwa_assembled")
#ribosome_dir <- file.path(data_dir, "ribosome_bwa_assembled")
output_dir <- "~/Desktop/pep_130618"

# Accession data
accessionData <- read.csv("/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/data/PeperomiaAccessionData.csv")
accessionData$label <- paste(accessionData$SampleID, accessionData$Genus, accessionData$Species, accessionData$Geography,sep = "_")

# Rename chloroplast files
chloroplast_files <- list.files(chloroplast_dir)
chloroplast_id <- str_split_fixed(chloroplast_files, pattern = "_", n = 2)[,1]

for(i in 1:length(chloroplast_files)){
  temp <- read.dna(file.path(chloroplast_dir, chloroplast_files[i]), "fasta")
  rownames(temp) <- accessionData$label[match(chloroplast_id[i], accessionData$SampleID)]
  write.dna(temp, file.path(output_dir, chloroplast_files[i]), "fasta")
}

ribosome_files <- list.files(ribosome_dir)
ribosome_id <- str_split_fixed(ribosome_files, pattern = "_", n = 2)[,1]

for(i in 1:length(chloroplast_files)){
  temp <- read.dna(file.path(ribosome_dir, ribosome_files[i]), "fasta")
  rownames(temp) <- accessionData$label[match(ribosome_id[i], accessionData$SampleID)]
  write.dna(temp, file.path(output_dir, ribosome_files[i]), "fasta")
}

head(accessionData)
peptaxa <- subset(accessionData, Genus == "Peperomia" & FinalLibrary.Y.N. == "Y")
sort(unique(as.vector(peptaxa$Species)))
