library(phytools)
library(ggtree)
phylo.dir <- "~/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/raxml_160618/"
data.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/data"

accessionData <- read.csv(file.path(data.dir, "PeperomiaAccessionData.csv"), stringsAsFactors = FALSE)
accessionData <- subset(accessionData, FinalLibrary.Y.N. == "Y")
accessionData$tiplabel <- paste(accessionData$SampleID,
                                accessionData$Genus,
                                accessionData$Species,
                                accessionData$Geography, sep = "_")
accessionData$tiplabel <- gsub(accessionData$tiplabel, pattern = " ", replacement = "_")
accessionData <- rbind(accessionData, rep(NA, ncol(accessionData)))
accessionData$tiplabel[nrow(accessionData)] <- "Piper_kadsura"
accessionData$Genus[nrow(accessionData)] <- "Piper"
accessionData$Species[nrow(accessionData)] <- "kadsura"

accessionData$newtiplabel <- paste(accessionData$Genus, accessionData$Species, accessionData$Geography)
rownames(accessionData) <- accessionData$label

accessionData$HigherGeography <- NA
#accessionData$HigherGeography[accessionData] Piper needs to be colorless
accessionData$HigherGeography[accessionData$Geography == "Australia" | accessionData$Geography == "New Zealand"] <- "Australia + New Zealand"
accessionData$HigherGeography[accessionData$Geography == "Africa" | accessionData$Geography == "Middle East"] <- "Africa"
accessionData$HigherGeography[accessionData$Geography == "Asia" | accessionData$Geography == "Japan"] <- "Asia"
accessionData$HigherGeography[accessionData$Geography == "Micronesia"] <- "Micronesia"
accessionData$HigherGeography[accessionData$Geography == "Hawaii"] <- "Hawaii"
accessionData$HigherGeography[accessionData$Geography %in% c("S Amer", "C Amer")] <- "S + C America"
accessionData$HigherGeography[accessionData$Geography %in% c("New Caledonia", "Fiji", "Samoa")] <- "W Pacific"
accessionData$HigherGeography[accessionData$Geography %in% c("Societies", "Australs", "Marquesas")] <- "S Pacific"

## PLOT PHYLOGENY

pepRAXML <- read.tree(file.path(phylo.dir, "RAxML_bipartitions.160618"))
pepRAXML_rooted <- phangorn::midpoint(pepRAXML)
#https://github.com/GuangchuangYu/ggtree/issues/89 (Describes some of the rooting issues I've been having in R)

# pepRAXML2 <- read.beast(file.path(phylo.dir, "bipartitions.nexus")) # file was opened in figtree, label was renamed as "bootstraps" and then saved as a nexus
# pepRAXML2@phylo <- phangorn::midpoint(pepRAXML2@phylo, node.lables = "label")

# pepRAXML2 <- read.raxml(file.path(phylo.dir, "RAxML_bipartitionsBranchLabels.160618"))
# pepRAXML2@data$bootstrap
# pepRAXML2@data$node
# max(pepRAXML2@phylo$edge)
# 
# pepRAXML2@phylo <- phangorn::midpoint(pepRAXML2@phylo)
# pepRAXML2@data$bootstrap
# pepRAXML2@data$node
# pepRAXML2@phylo$edge

# testPhylo_support <- ggtree(pepRAXML_rooted) + geom_text2(aes(label=label), hjust=-.2)
# read.tree import node labels as text, for some reason read.raxml assumes they are node labels and not the bootstraps
# only show supports > 80?

# pepRAXML2@phylo$tip.label # 122 tips
# pepRAXML2@phylo$tip.label %in% accessionData$tiplabel # 121 tips

library(ggplot2)
testPhylo_support <- ggtree(pepRAXML_rooted, size = 0.2) %<+%
  accessionData[c("tiplabel","newtiplabel", "HigherGeography")] +
  #geom_text2(aes(label=label, subset=!is.na(as.numeric(label)) & as.numeric(label)>80 ), size = 3, hjust = -.2) + # plots bootstraps
  geom_nodepoint(aes(subset=!is.na(as.numeric(label)) & as.numeric(label)>80), color = "grey30", alpha = 0.5, size= 3) + # adds a dot for nodes > 80 bootstrap
  geom_tiplab(aes(label=newtiplabel, color = HigherGeography), size = 1.5) +
  xlim(0, 0.045) +
  scale_color_manual(values = colBiogeog) +
  geom_treescale(x = 0.025, y = -1, offset = -2)

sort(unique(accessionData$HigherGeography))
col1 <- rgb(120, 135, 189, maxColorValue = 255)
col2 <- rgb(233, 166, 63, maxColorValue = 255)
col3 <- rgb(142, 36, 118, maxColorValue = 255)
col4 <- rgb(133, 196, 211, maxColorValue = 255)
col5 <- rgb(55, 125, 155, maxColorValue = 255)
col6 <- rgb(171, 49, 44, maxColorValue = 255)
col7 <- rgb(121, 176, 84, maxColorValue = 255)
col8 <- rgb(241, 216, 110, maxColorValue = 255)
colBiogeog <- c(col4, col7, col2, col1, col5, col6, col3, col8)
plot(1:8, col = c(col1, col2, col3, col4, col5, col6, col7, col8), cex = 3, pch = 16)
  
ggsave(testPhylo_support, width = 8, height = 12, filename = "~/Desktop/pepSupport.pdf")


geom_tiplab(aes(color = Geography), size =1.5) +

# For this method to work, the FIRST column must be matched with the tip labels
  
  
testPhylo_node <- ggtree(pepRAXML_rooted) + geom_tiplab(size = 1.2) + geom_text(aes(label=node))





ggsave(testPhylo_node, width = 20, height = 20, filename = "~/Desktop/pepNode.pdf")
