## Packages
library(phytools)

## Directories
main.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/"
data.dir <- file.path(main.dir, "data")
tree.dir <- file.path(main.dir, "output")

## Import accession data
test <- read.csv(file.path(data.dir, "PeperomiaAccessionData.csv"))


pep_tree <- read.newick(file.path(tree.dir, "PS_gene_map.tre"))
pep_tree
