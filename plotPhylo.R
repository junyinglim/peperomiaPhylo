## Plot phylogenies

## Packages ============
library(phytools); library(tibble)
library(ggtree); library(tidytree)
library(ggplot2); library(gridExtra); library(grid)
library(stringr); library(deeptime)

## Directories ============
phylo.dir <- "~/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/2018-09-04/iqtree"
data.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/data"
main.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/"
output.dir <- file.path(main.dir, "output")
fig.dir <- file.path(main.dir, "figures")

## Import accessionData ============
accessionData <- read.csv(file.path(data.dir, "PeperomiaAccessionData.csv"), stringsAsFactors = FALSE)
accessionData <- subset(accessionData, FinalLibrary.Y.N. == "Y")

accessionData$tiplabel <- accessionData$SampleID
# accessionData$tiplabel <- paste(accessionData$SampleID,
#                                 accessionData$Genus,
#                                 accessionData$Species,
#                                 accessionData$Geography, sep = "_")

accessionData <- rbind(accessionData, rep(NA, ncol(accessionData)))
accessionData$tiplabel[nrow(accessionData)] <- "Piper_kadsura"
accessionData$Genus[nrow(accessionData)] <- "Piper"
accessionData$Species[nrow(accessionData)] <- "kadsura"

# Small fixes
accessionData$Genus[accessionData$tiplabel == "PEZ-247"]  <- "Piper"
accessionData$newtiplabel <- paste(accessionData$Genus, accessionData$Species, accessionData$Geography)
rownames(accessionData) <- accessionData$label

accessionData$newtiplabel[accessionData$newtiplabel == "Piper kadsura NA"] <- "Piper kadsura"
accessionData$newtiplabel[accessionData$tiplabel == "PEZ-180_Peperomia_quadrifolia_C_Amer"] <- "Peperomia sp. C Amer"

# Define higher geography categories
accessionData$HigherGeography <- NA
accessionData$HigherGeography[accessionData$Geography == "Australia" | accessionData$Geography == "New Zealand"] <- "Australia + New Zealand"
accessionData$HigherGeography[accessionData$Geography == "Africa" | accessionData$Geography == "Middle East"] <- "Africa"
accessionData$HigherGeography[accessionData$Geography == "Asia" | accessionData$Geography == "Japan"] <- "Asia"
accessionData$HigherGeography[accessionData$Geography == "Micronesia"] <- "Micronesia"
accessionData$HigherGeography[accessionData$Geography == "Hawaii"] <- "Hawaii"
accessionData$HigherGeography[accessionData$Geography %in% c("Juan Fernandez","S Amer", "C Amer")] <- "S + C America"
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

## Summarize dated phylogenies ============
treePLtrees <- list.files(output.dir, pattern = "treePL_[0-9]{1,3}.out$")
treePLtrees_brtimes <- lapply(treePLtrees, FUN = function(x) { tree <- read.tree(file.path(output.dir, x)); return(branching.times(tree))})
treePLtrees_brtimes2 <- do.call("cbind", treePLtrees_brtimes)
treePLtrees_brtimes_range <- apply(X = treePLtrees_brtimes2, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975))

treePLtrees_mean <- read.tree(file.path(output.dir, "treePL_mean.out"))
treePLtrees_brtimes_mean <- branching.times(treePLtrees_mean)

brtimesList <- list()
for(i in 1:ncol(treePLtrees_brtimes_range)){
  brtimesList[[i]] <- c(treePLtrees_brtimes_range[1,i], treePLtrees_brtimes_range[2,i])
}

brtimesTipEmpty <- vector("list",length(treePLtrees_mean$tip.label))
brtimesTipEmpty[1:length(treePLtrees_mean$tip.label)] <- NA
z <- append(brtimesTipEmpty, brtimesList)
treePLtrees_mean2 <- as_data_frame(treePLtrees_mean)
treePLtrees_mean3 <- as.treedata(treePLtrees_mean2)
treePLtrees_mean3@data$"height_range"=z
treePLtrees_mean3@data$"height"=c(rep(NA,120), treePLtrees_brtimes_mean)

treePLtrees_mean3@phylo$node.label <- (Ntip(treePLtrees_mean3)+1):(Ntip(treePLtrees_mean3)+Nnode(treePLtrees_mean3))


test <- revts(ggtree(treePLtrees_mean3) +
  coord_cartesian(xlim = c(-80,0), ylim = c(-10, Ntip(treePLtrees_mean3)+2), expand = FALSE) +
  geom_range(range = 'height_range', branch.length = "height", col = "navyblue", alpha = 0.5, size = 2) + #geom_nodelab() +
#  geom_cladelabel(node = 130, label = "Hawaiian endemics A", offset = -1, angle = 90, ) +
  scale_x_continuous(breaks=seq(-80,0, 10), labels=abs(seq(-80,0,10)),name = "Millions of years (Ma)") +
  geom_vline(aes(xintercept = -66), alpha = 0.5, col = "grey30", linetype = "dashed") +
  geom_vline(aes(xintercept = -56), alpha = 0.5, col = "grey30", linetype = "dashed") +
  geom_vline(aes(xintercept = -33.9), alpha = 0.5, col = "grey30", linetype = "dashed") +
  geom_vline(aes(xintercept = -23.03), alpha = 0.5, col = "grey30", linetype = "dashed") +
  geom_vline(aes(xintercept = -5.333), alpha = 0.5, col = "grey30", linetype = "dashed") +
  geom_vline(aes(xintercept = -2.58), alpha = 0.5, col = "grey30", linetype = "dashed") +
  theme_tree2() )
test2 <- gggeo_scale(test, dat = "epochs", neg = TRUE, abbrv = TRUE, size = 2, lab = FALSE)

ggsave(test2, filename = file.path(fig.dir, "datedPepPhy.pdf"), device = cairo_pdf, width = 5, height = 5)

# ggtree(treePLtrees_mean) + geom_segment2(aes(subset=!isTip, xend = c(rep(NA, 120), treePLtrees_brtimes_range[1,]), yend = c(rep(NA, 120), treePLtrees_brtimes_range[2,]))) + scale_x_reverse(breaks=seq(-50, 0, by = 10), labels=abs(seq(-50,0,by = 10)),expand = c(0,0),name = "Millions of years (Ma)")

## Plot phylogeny ==============
# Import phylogeny -----------
pepML <- read.nexus(file.path(phylo.dir, "output.treefile_rooted"))
#pepRAXML_rooted <- phangorn::midpoint(pepRAXML)
#https://github.com/GuangchuangYu/ggtree/issues/89 (Describes some of the rooting issues I've been having in R)
pepML$tip.label <- gsub(pepML$tip.label, pattern = "_bwa", replacement = "")

# Correct one of the mislabels
pepML$tip.label[pepML$tip.label == "PEZ-294"] <- "filler" #"PEZ-298_Peperomia_adamsonii_Marquesas"
pepML$tip.label[pepML$tip.label == "PEZ-298"] <- "PEZ-294"
pepML$tip.label[pepML$tip.label == "filler"] <- "PEZ-298" #"PEZ-298_Peperomia_adamsonii_Marquesas"

# Tidy up bootstrap labels
pepML$node.label <- gsub(pepML$node.label, pattern = "\"|\"", replacement = "")
pepML$node.label <- str_split_fixed(pepML$node.label, pattern = "/", n = 2)[,2] #bootstrap

# Drop a labelling error
pepML <- drop.tip(pepML, tip = "PEZ-211") # Peperomia membranaceae

# Define colors-----------
#sort(unique(accessionData$HigherGeography))
col1 <- rgb(120, 135, 189, maxColorValue = 255)
col2 <- rgb(233, 166, 63, maxColorValue = 255)
col3 <- rgb(142, 36, 118, maxColorValue = 255)
col4 <- rgb(133, 196, 211, maxColorValue = 255)
col5 <- rgb(55, 125, 155, maxColorValue = 255)
col6 <- rgb(171, 49, 44, maxColorValue = 255)
col7 <- rgb(121, 176, 84, maxColorValue = 255)
col8 <- rgb(241, 216, 110, maxColorValue = 255)
colBiogeog <- c(col4, col7, col2, col1, col5, col8, col6, col3, "grey60")
#plot(1:8, col = c(col1, col2, col3, col4, col5, col6, col7, col8), cex = 3, pch = 16)

# Create legend -----------
legendDF <- data.frame(cat=sort(unique(accessionData$HigherGeography)), x = 1:9, y = 1.9)
testLegend <- ggplot(data = legendDF) + geom_point(aes(y = y, x = x, color = cat), size = 3) + scale_colour_manual(values = colBiogeog) + theme(legend.key = element_blank(), legend.title = element_blank())

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

plotLegend <- g_legend(testLegend)
ggsave(plotLegend, filename = file.path(output.dir, "biogeogLegend.pdf"), device = cairo_pdf, width = 2, height = 4)

# Plot phylogeny -----------
# For this method to work, the FIRST column must be matched with the tip labels
phyloSupport <- ggtree(pepML, size = 0.2, ladderize = TRUE) %<+%
  accessionData[c("tiplabel","newtiplabel", "HigherGeography")] +
  #geom_nodelab(aes(subset=!is.na(as.numeric(label)) & as.numeric(label)<80), size = 1.5) +
  geom_nodepoint(aes(subset=!is.na(as.numeric(label)) & as.numeric(label)>80), color = "grey30", alpha = 0.5, size= 1.5) + # adds a dot for nodes > 80 bootstrap
  #geom_hilight(node = 214, fill = "grey90", extendto = 0.035) +
  #geom_hilight(node = 132, fill = "grey90", extendto = 0.035) +
  geom_tiplab(aes(label=newtiplabel, color = HigherGeography), size = 1.5) +
  xlim(0, 0.07) +
  #geom_cladehi(node = 214, label = "Hawaiian endemics A") +
  #geom_cladehi(node = 132, label = "Hawaiian endemics B") +
  scale_color_manual(values = colBiogeog) +
  #geom_treescale(x = 0, y = 0, offset = -2) +
  theme(panel.background = element_rect(fill = "transparent"))

#phyloNode <- ggtree(pepRAXML_rooted, size = 0.2, ladderize = TRUE) + 
#  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 1)

#pepPhylogram <- arrangeGrob(phyloSupport, g_legend(testLegend), nrow = 1, ncol = 2, widths = c(8, 2)) 
ggsave(phyloSupport, width = 8, height = 12, filename = file.path(fig.dir, "phylo_support.pdf"), device = cairo_pdf)

## PLOT THE DATED PHYLOGENY =======
datedPhy <- read.tree(file.path(phylo.dir, "dated.tre"))

# Correct some labelling errors
datedPhy$tip.label[datedPhy$tip.label == "PEZ-294_Peperomia_tetraphylla_Australia"] <- "filler"
datedPhy$tip.label[datedPhy$tip.label == "PEZ-298_Peperomia_adamsonii_Marquesas"] <- "PEZ-294_Peperomia_tetraphylla_Australia"
datedPhy$tip.label[datedPhy$tip.label == "filler"] <- "PEZ-298_Peperomia_adamsonii_Marquesas"

datedPhy <- drop.tip(datedPhy, tip = "PEZ-211_Peperomia_membranacea_Hawaii")

datedPhy <- rotateConstr(datedPhy, pepRAXML_rooted$tip.label) # ensure that tips are in the same order

datedPepPhy <- drop.tip(datedPhy, tip = c("PEZ-204_Piper_ponapense_Micronesia", "PEZ-247_Peperomia_sp._New_Caledonia", "Piper_kadsura", "PEZ-164_Macropiper_puberulum_Samoa")) #remove piper from tree

datedPlot <- ggtree(datedPepPhy, size = 0.2) +
  theme_tree2()

# datedPlot2 <- revts(datedPlot) + scale_x_continuous(breaks=seq(-50, 0, by = 10),
#                                                    labels=abs(seq(-50,0, by = 10)),
#                                                    limits = c(-55,0),
#                                                    expand = c(0,0) 
#                                                   )

datedPlot2 <- revts(datedPlot) + scale_x_reverse(breaks=seq(-50, 0, by = 10),
                                                labels=abs(seq(-50,0,by = 10)),
                                                expand = c(0,0),
                                                name = "Millions of years (Ma)") +
  geom_vline(xintercept = seq(-50,0, by = 5), size = 0.2, alpha = 0.7)


                                
datedPlot3 <- ggtree(datedPhy, size = 0.2) + scale_x_reverse()

# cowplot::plot_grid(testPhylo_support, datedPlot3, ncol=2)

ggsave(datedPlot2, filename = file.path(phylo.dir, "datedPlot.pdf"), height = 12, width = 6)
