library(phytools)

main.dir <- "~/Dropbox/Projects/2015/Peperomia/data/output/"
data.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/data"

accessionData <- read.csv(file.path(data.dir, "PeperomiaAccessionData.csv"))
accessionData$tiplabel <- paste(accessionData$SampleID, accessionData$Species, accessionData$Geography, sep = "_")

tree <- read.nexus(file.path(main.dir, "PS_gene_map.tre"))

tree$tip.label <- accessionData$tiplabel[match(tree$tip.label, accessionData$SampleID)]

tree$tip.label[33] <- "Piper_kadsura"

write.tree(tree, file.path(main.dir, "pepPhy.tre"))
