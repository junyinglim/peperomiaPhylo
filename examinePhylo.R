## PACKAGES ============
library(phytools)
library(ggtree)
library(ggplot2)
library(gridExtra)
library(grid)

## DIRECTORIES ============
data.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/data"
main.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/"
phylo.dir <- "~/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/2018-09-04/iqtree/"
# IQ-TREE run fully partitioned, all 165 protein-coding and non-coding regions

phylo.dir <- "~/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/2018-09-04/raxml/"

## IMPORT FILES ============
# Import tree
pepIQTREE_rooted <- read.nexus(file.path(phylo.dir, "output.treefile_rooted"))
# NOTE: original file needs to be opened in figtree, rooted manually, saved as a nexus file, and the square brackets removed by regular expr.
# WHY THE TROUBLE? midpoint function makes mistakes when renaming nodes
pepIQTREE_rooted <- read.tree(file.path(phylo.dir, "RAxML_bipartitions.raxml_2018-09-04"))

pepIQTREE_rooted$tip.label <- gsub(pepIQTREE_rooted$tip.label, pattern = "_bwa", replacement = "")

# Import accession data
accessionData <- read.csv(file.path(data.dir, "PeperomiaAccessionData.csv"), stringsAsFactors = FALSE)
accessionData <- subset(accessionData, FinalLibrary.Y.N. == "Y")
accessionData$tiplabel <- accessionData$SampleID

accessionData <- rbind(accessionData, rep(NA, ncol(accessionData)))
accessionData$tiplabel[nrow(accessionData)] <- "Piper_kadsura"
accessionData$Genus[nrow(accessionData)] <- "Piper"
accessionData$Species[nrow(accessionData)] <- "kadsura"

accessionData$Genus[accessionData$tiplabel == "PEZ-247"]  <- "Piper" # Small fix to genus determination

# Create new, tip labels for display
accessionData$newtiplabel <- paste(accessionData$Genus, accessionData$Species, accessionData$Geography, accessionData$SampleID)
rownames(accessionData) <- accessionData$label

# Small fixes
accessionData$newtiplabel[accessionData$newtiplabel == "Piper kadsura NA NA"] <- "Piper kadsura"
accessionData$newtiplabel[accessionData$tiplabel == "PEZ-180"] <- "Peperomia sp. C Amer PEZ-180"

# Define higher geography categories
accessionData$HigherGeography <- NA
accessionData$HigherGeography[accessionData$Geography == "Australia" | accessionData$Geography == "New Zealand"] <- "Australia + New Zealand"
accessionData$HigherGeography[accessionData$Geography == "Africa" | accessionData$Geography == "Middle East"] <- "Africa"
accessionData$HigherGeography[accessionData$Geography == "Asia" | accessionData$Geography == "Japan"] <- "Asia"
accessionData$HigherGeography[accessionData$Geography == "Micronesia"] <- "Micronesia"
accessionData$HigherGeography[accessionData$Geography == "Hawaii"] <- "Hawaii"
accessionData$HigherGeography[accessionData$Geography %in% c("S Amer", "C Amer", "Juan Fernandez")] <- "S + C America"
accessionData$HigherGeography[accessionData$Geography %in% c("New Caledonia", "Fiji", "Samoa")] <- "New Caledonia + Fiji + Samoa"
accessionData$HigherGeography[accessionData$Geography %in% c("Pitcairns", "Societies", "Australs", "Marquesas")] <- "S Pacific"
accessionData$HigherGeography[accessionData$Genus %in% c("Piper", "Macropiper")] <- "Outgroups"

accessionData$HigherGeography <- factor(accessionData$HigherGeography,
                                        levels = c("Africa",
                                                   "Asia",
                                                   "Australia + New Zealand",
                                                   "Hawaii",
                                                   "Micronesia",
                                                   "New Caledonia + Fiji + Samoa",
                                                   "S + C America",
                                                   "S Pacific",
                                                   "Outgroups"))


# Correct one of the mislabels
pepIQTREE_rooted$tip.label[pepIQTREE_rooted$tip.label == "PEZ-294"] <- "filler" #"PEZ-298_Peperomia_adamsonii_Marquesas"
pepIQTREE_rooted$tip.label[pepIQTREE_rooted$tip.label == "PEZ-298"] <- "PEZ-294"
pepIQTREE_rooted$tip.label[pepIQTREE_rooted$tip.label == "filler"] <- "PEZ-298" #"PEZ-298_Peperomia_adamsonii_Marquesas"

# Drop a labelling error (perhaps leave as is?)
#pepIQTREE_rooted <- drop.tip(pepIQTREE_rooted, tip = "PEZ-211_Peperomia_membranacea_Hawaii")

## PLOT PHYLOGENY ===============

# Define biogeographic colors
sort(unique(accessionData$HigherGeography))
col1 <- rgb(120, 135, 189, maxColorValue = 255)
col2 <- rgb(233, 166, 63, maxColorValue = 255)
col3 <- rgb(142, 36, 118, maxColorValue = 255)
col4 <- rgb(133, 196, 211, maxColorValue = 255)
col5 <- rgb(55, 125, 155, maxColorValue = 255)
col6 <- rgb(171, 49, 44, maxColorValue = 255)
col7 <- rgb(121, 176, 84, maxColorValue = 255)
col8 <- rgb(241, 216, 110, maxColorValue = 255)
colBiogeog <- c(col4, col7, col2, col1, col5, col8, col6, col3, "grey60")

# Clean up node labels
pepIQTREE_rooted$node.label <- gsub(pepIQTREE_rooted$node.label, pattern = "\"|\"", replacement = "")

phyloSupport <- ggtree(pepIQTREE_rooted, size = 0.2, ladderize = TRUE) %<+%
  accessionData[c("tiplabel","newtiplabel", "HigherGeography")] +
  geom_nodelab(size = 1) +
  geom_tiplab(aes(label=newtiplabel, color = HigherGeography), size = 1.5) +
  #xlim(0, 0.045) +
  scale_color_manual(values = colBiogeog) +
  geom_treescale(x = 0, y = 0, offset = -2) +
  theme(panel.background = element_rect(fill = "transparent"))

ggsave(phyloSupport, width = 18, height = 12, filename = file.path(phylo.dir, "pepPhylo_support.pdf"), device = cairo_pdf)



pepIQTREE_rooted_trunc <- drop.tip(pepIQTREE_rooted, tip = c("PEZ-247", "PEZ-164", "PEZ-204", "Piper_kadsura"))
phyloSupportTrunc <- ggtree(pepIQTREE_rooted_trunc, size = 0.2, ladderize = TRUE) %<+%
  accessionData[c("tiplabel","newtiplabel", "HigherGeography")] +
  geom_nodelab(size = 1) +
  geom_tiplab(aes(label=newtiplabel, color = HigherGeography), size = 1.5) +
  #xlim(0, 0.045) +
  scale_color_manual(values = colBiogeog) +
  geom_treescale(x = 0, y = 0, offset = -2) +
  theme(panel.background = element_rect(fill = "transparent"))
ggsave(phyloSupportTrunc, width = 18, height = 12, filename = file.path(phylo.dir, "pepPhylo_supportTrunc.pdf"), device = cairo_pdf)
