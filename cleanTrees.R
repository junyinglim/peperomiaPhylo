## Packages
library(phytools)

## Directories
main.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/"
data.dir <- file.path(main.dir, "data")
tree.dir <- file.path(main.dir, "output")

## Import accession data
accessionData <- read.csv(file.path(data.dir, "PeperomiaAccessionData.csv"))
accessionData$tiplabel <- paste(accessionData$SampleID, accessionData$Species, accessionData$Geography, sep = "_")

pep_tree <- read.newick(file.path(tree.dir, "PS_gene_map.tre"))
pep_tree



pepPhy<- read.tree("~/Desktop/peperomiaRAXML/RAxML_bipartitions.pep_nopart")
index <- match(pepPhy$tip.label, test$SampleID)

pepPhy$tip.label <- accessionData$tiplabel[match(pepPhy$tip.label, accessionData$SampleID)]
pepPhy$tip.label[which(is.na(index))] <- "Piper_kadsura"
write.tree(pepPhy, "~/Desktop/peperomiaRAXML/pepPhy.tre")


