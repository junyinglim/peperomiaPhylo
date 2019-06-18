## Plot phylogenies

## Packages ============
library(phytools); library(tibble)
library(ggtree); library(tidytree)
library(ggplot2); library(gridExtra); library(grid)
library(stringr); library(deeptime)
library(RColorBrewer); library(reshape2)

## Directories ============
phylo.dir <- "~/Dropbox/Peperomia/Pacific/genomeskimming_data/pepPhyloRuns/2018-09-04/iqtree"
main.dir <- "~/Dropbox/Peperomia/Pacific/peperomiaPhylo/"
data.dir <- file.path(main.dir, "data")
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
accessionData$Geography[nrow(accessionData)] <- ""

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
  geom_range(range = 'height_range', branch.length = "height", col = "navyblue", alpha = 0.5, size = 2) +
  geom_nodelab(size = 2) +
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

treePLtrees_mean3@data$height_range[130] # crown-age radiation A
treePLtrees_mean3@data$height[130]

treePLtrees_mean3@data$height_range[209] # crown-age radiation B
treePLtrees_mean3@data$height[209]

treePLtrees_mean3@data$height_range[129] # stem-age radiation A
treePLtrees_mean3@data$height[129]

treePLtrees_mean3@data$height_range[208] # stem-age radiation B
treePLtrees_mean3@data$height[208]

## Summarize fossilized birth-death trees
#mapTree <- read.beast(file.path(output.dir, "ucln_map.tree"))
mccTree <- read.beast(file.path(output.dir, "ucln_mcc.tree"))
mccTree_extant <- treeio::drop.tip(mccTree, tip = c("Saururus_aquilae", "Saururus_tuckerae", "Saururus_bilobatus", "Saururus_stoobensis", "Houttuynia_bavarica", "Saururopsis_niponensis", "Hexagyne_philippiana", "Lactoripollenites_africanus", "Aristolochia_austriaca", "Piper_margaritae", "Piper_bartlingianum"))
mccTree_extant@data$age <- c(rep(NA, Ntip(mccTree_extant)),branching.times(mccTree_extant@phylo))

mccTreePlot <- revts(ggtree(mccTree_extant) + geom_tiplab(size = 2) + 
                       coord_cartesian(xlim = c(-140,40), ylim = c(-10, Ntip(mccTree_extant)+2), expand = FALSE) + 
        geom_range("age_0.95_HPD", branch.length = "age", color = "red", size = 1.8, alpha = 0.5) + #geom_nodelab() +
      scale_x_continuous(breaks=seq(-140,0, 10), labels=abs(seq(-140,0,10)),name = "Millions of years (Ma)") +
  theme_tree2() )

mccTreePlot2 <- gggeo_scale(mccTreePlot, dat = "periods", neg = TRUE, abbrv = FALSE, size = 2)
ggsave(mccTreePlot2, filename = file.path(fig.dir, "datedPiperales.pdf"), device = cairo_pdf, width = 7, height = 7)

FBDlogFiles <- list.files(output.dir, pattern = "ucln_run_[0-9]{1,2}\\.log")
FBDlogDF <- lapply(file.path(output.dir, FBDlogFiles), FUN = read.table, header = TRUE)
FBDlogDF_burnin <- lapply(FBDlogDF, FUN = function(x){x[round(0.25*nrow(x)):nrow(x),]})
FBDlogDF_total <- do.call("rbind",FBDlogDF_burnin)
quantile(FBDlogDF_total$crownAge_pippep, probs = c(0.025, 0.975))
mean(FBDlogDF_total$crownAge_pippep)

mccTree@data$age_0.95_HPD[mccTree@data$node == 122]
branching.times(mccTree@phylo)["122"]

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

# Plot phylogeny -----------
# Define colours for biogeographic regions
berkcol3 = c("#003262", "#3B7EA1", "#9BBEA9", "#00B0DA", "#00A598", brewer.pal(9, "BuGn")[8], "#CFDD45","#859438", "#FDB515", brewer.pal(9, "YlOrRd")[5], "#ED4E33", "#C4820E", "#D9661F","#6C3302","grey80")
#plot(1:15, col =berkcol3, pch = 16, cex=3)

# For this method to work, the FIRST column must be matched with the tip labels
accessionData$label <- accessionData$SampleID
accessionData$binom <- paste(accessionData$Genus, accessionData$Species)
accessionData$binom <- paste0("`", accessionData$binom, "`")
rownames(accessionData) <- accessionData$tiplabel
accessionData_subset <- subset(accessionData, tiplabel %in% pepML$tip.label)
accessionData_subset <- accessionData_subset[match(accessionData_subset$tiplabel, pepML$tip.label),]

accessionData_subset$HigherGeography <- factor(accessionData_subset$HigherGeography,
                                               levels = c("Hawaii",
                                                          "S Pacific",
                                                          "New Caledonia + Fiji + Samoa",
                                                          "Asia",
                                                          "Australia + New Zealand",
                                                          "Micronesia",
                                                          "Africa",
                                                          "S + C America",
                                                          "Outgroups"))
accessionData_subset2 <- dcast(data = accessionData_subset, formula = tiplabel~HigherGeography, value.var = "HigherGeography")
rownames(accessionData_subset2) <- accessionData_subset2$tiplabel
target_col <- levels(accessionData_subset$HigherGeography)

# Plot phylogeny with tip labels
phyloSupport <- ggtree(pepML, size = 0.2, ladderize = TRUE) %<+%
  accessionData_subset[c("label","tiplabel","binom","HigherGeography")] +
  geom_nodepoint(aes(subset=!is.na(as.numeric(label)) & as.numeric(label)>80 & as.numeric(label)<=90), color = "black", fill = "grey80", size= 1.5, shape = 21) + 
  geom_nodepoint(aes(subset=!is.na(as.numeric(label)) & as.numeric(label)>90 & as.numeric(label)<=100), color = "black", fill = "white", size= 1.5, shape = 21) +
  geom_tiplab(aes(label=paste0("italic(",binom,")")), size = 1.5, parse = TRUE) +
  #xlim(0, 0.07) +
  #scale_color_manual(values = colBiogeog) +
  ylim(-2,119) +
  geom_treescale(x = 0, y = 0, offset = -1.8) +
  theme(panel.background = element_rect(fill = "transparent"))
ggsave(phyloSupport, width = 8, height = 12, filename = file.path(fig.dir, "phylo_support.pdf"))

phyloBiogeog <- gheatmap(p = phyloSupport,  data = accessionData_subset2[target_col], width = 0.08, color = "transparent", offset= 0.007, colnames = F) +
  scale_fill_manual(breaks = target_col,
                    values = c(berkcol3[c(1,3,4,6,7,8,15,9,11)]),
                    labels = target_col) +
  theme(legend.position = "right", legend.title = element_blank())
ggsave(phyloBiogeog + theme(legend.position = "none"),
       width = 9, height = 12, filename = file.path(fig.dir, "phylo_biogeog.pdf"))


phyloBiogeog_legend <- cowplot::get_legend(phyloBiogeog + theme(legend.spacing.x = unit(0.3,"cm"), legend.background = element_blank()))
ggsave(phyloBiogeog_legend, filename = file.path(fig.dir, "biogeogLegend.pdf"), device = cairo_pdf, width = 3, height = 4, bg = "transparent")
#phyloNode <- ggtree(pepRAXML_rooted, size = 0.2, ladderize = TRUE) + 
#  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 1)

## QUARTET SAMPLING RESULTS =======
qs_qc <- read.newick("~/Dropbox/Peperomia/Pacific/genomeskimming_data/pepPhyloRuns/2018-09-04/qs/RESULT.labeled.tre.qc.rooted.nex")
qs_qi <- read.newick("~/Dropbox/Peperomia/Pacific/genomeskimming_data/pepPhyloRuns/2018-09-04/qs/RESULT.labeled.tre.qi.rooted.nex")

# Correct one of the mislabels
cleanUpTreeTips <- function(tree){
  #tree <- phangorn::midpoint(tree, node.labels= "label")
  tree$tip.label[tree$tip.label == "PEZ-294"] <- "filler" #"PEZ-298_Peperomia_adamsonii_Marquesas"
  tree$tip.label[tree$tip.label == "PEZ-298"] <- "PEZ-294"
  tree$tip.label[tree$tip.label == "filler"] <- "PEZ-298" #"PEZ-298_Peperomia_adamsonii_Marquesas"
  tree <- drop.tip(tree, tip = "PEZ-211") # Peperomia membranaceae
  return(tree)
}

# Focus on qi and qc
qs_qi_clean <- cleanUpTreeTips(qs_qi)
qs_qc_clean <- cleanUpTreeTips(qs_qc)
qs_qi_clean$node.label <- gsub(qs_qi_clean$node.label, pattern = "qi=", replacement = "")
qs_qc_clean$node.label <- gsub(qs_qc_clean$node.label, pattern = "qc=", replacement = "")

library(viridis)
qs_qi_plot <- ggtree(qs_qi_clean)  %<+%
  accessionData[c("tiplabel","binom","HigherGeography")] +
  geom_tiplab(aes(label=paste0("italic(",binom,")")), size = 1.5, parse = TRUE) +
  geom_nodepoint(aes(colour = as.numeric(label))) +
  xlim(0, 0.06) +
  scale_colour_viridis(name = "Quartet Informativeness") +
  theme(legend.position = "bottom")

qs_qc_plot <- ggtree(qs_qc_clean)  %<+%
  accessionData[c("label","tiplabel","binom","HigherGeography")] +
  geom_tiplab(aes(label=paste0("italic(",binom,")")), size = 1.5, parse = TRUE) +
  geom_nodepoint(aes(colour = as.numeric(label))) +
  xlim(0, 0.06) +
  scale_colour_viridis(name = "Quartet Concordance") +
  theme(legend.position = "bottom")

ggsave(plot = qs_qi_plot, file.path(fig.dir, "pepML_qs_qi.pdf"), width = 8, height = 12)
ggsave(plot = qs_qc_plot, file.path(fig.dir, "pepML_qs_qc.pdf"), width = 8, height = 12)

qs_plots <- plot_grid(qs_qi_plot, qs_qc_plot, labels = c("(a)", "(b)"), nrow = 2, label_size = 20)

ggsave(qs_plots, filename = file.path(fig.dir, "pepML_qs_comb.pdf"), width = 8, height = 12)

# QC = Quartet concordance - how often quartet is preferred over other config (1 = concord., 0 = equivocal, <0 = disc. > conc.)
# QD = Quartet differential - are discordant frequencies equal or skewed? (1 = equal, 0 = one discordant is completely dominant)
# QI = Quartet informativeness - what prop. of replicates exceeded likelihood differential (1 - all infomrative, 0 = none informative)
# QF = Quartet fidelity - when a taxon is sampled, how often does it produce a concordant topology (1 = all cncordant, 0.1 = 10% concordant, 0 = none concordant) [ used to identify misidentification / ingroup errors ]

## SUMMARY STATISTICS =======
z <- subset(accessionData, tiplabel %in% pepML$tip.label) # 118 accessions (incl 3 outgroups)
table(subset(accessionData, tiplabel %in% pepML$tip.label)$HigherGeography) # 97 accessions

# 2 Piper
length(pepML$tip.label)
z2 <- subset(z, ! HigherGeography %in% c("Africa", "Asia", "Outgroups", "S + C America"))
sort(unique(paste0(z2$Genus, z2$Species)))


readDepthSummary <- read.csv(file.path(main.dir, "readDepthSummary.csv"))
readDepthSummary_subset <- subset(readDepthSummary, sample_ID %in% pepML$tip.label)
range(readDepthSummary_subset$meanCoverage)

## 
