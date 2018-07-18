
## PACKAGES ============
library(phytools)
library(ggtree)
library(ggplot2)
library(gridExtra)
library(grid)

## DIRECTORIES ============
phylo.dir <- "~/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/raxml_160618/"
data.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/data"
main.dir <- "~/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/"

# Import accessionData
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

# NOTE: SMALL FIX TO GENUS DET.
accessionData$Genus[accessionData$tiplabel == "PEZ-247_Peperomia_sp._New_Caledonia"]  <- "Piper"

accessionData$newtiplabel <- paste(accessionData$Genus, accessionData$Species, accessionData$Geography)
rownames(accessionData) <- accessionData$label

# Another small fix to remove biogeography 
accessionData$newtiplabel[accessionData$newtiplabel == "Piper kadsura NA"] <- "Piper kadsura"
accessionData$newtiplabel[accessionData$tiplabel == "PEZ-180_Peperomia_quadrifolia_C_Amer"] <- "Peperomia sp. C Amer"

# Define higher geography categories
accessionData$HigherGeography <- NA
accessionData$HigherGeography[accessionData$Geography == "Australia" | accessionData$Geography == "New Zealand"] <- "Australia + New Zealand"
accessionData$HigherGeography[accessionData$Geography == "Africa" | accessionData$Geography == "Middle East"] <- "Africa"
accessionData$HigherGeography[accessionData$Geography == "Asia" | accessionData$Geography == "Japan"] <- "Asia"
accessionData$HigherGeography[accessionData$Geography == "Micronesia"] <- "Micronesia"
accessionData$HigherGeography[accessionData$Geography == "Hawaii"] <- "Hawaii"
accessionData$HigherGeography[accessionData$Geography %in% c("S Amer", "C Amer")] <- "S + C America"
accessionData$HigherGeography[accessionData$Geography %in% c("New Caledonia", "Fiji", "Samoa")] <- "New Caledonia + Fiji + Samoa"
accessionData$HigherGeography[accessionData$Geography %in% c("Societies", "Australs", "Marquesas")] <- "S Pacific"
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

## PLOT PHYLOGENY ==============
pepRAXML <- read.tree(file.path(phylo.dir, "RAxML_bipartitions.160618"))
pepRAXML_rooted <- phangorn::midpoint(pepRAXML)
#https://github.com/GuangchuangYu/ggtree/issues/89 (Describes some of the rooting issues I've been having in R)

# Correct one of the mislabels
pepRAXML_rooted$tip.label[pepRAXML_rooted$tip.label == "PEZ-294_Peperomia_tetraphylla_Australia"] <- "filler" #"PEZ-298_Peperomia_adamsonii_Marquesas"
pepRAXML_rooted$tip.label[pepRAXML_rooted$tip.label == "PEZ-298_Peperomia_adamsonii_Marquesas"] <- "PEZ-294_Peperomia_tetraphylla_Australia"
pepRAXML_rooted$tip.label[pepRAXML_rooted$tip.label == "filler"] <- "PEZ-298_Peperomia_adamsonii_Marquesas" #"PEZ-298_Peperomia_adamsonii_Marquesas"

# Drop a labelling error
pepRAXML_rooted <- drop.tip(pepRAXML_rooted, tip = "PEZ-211_Peperomia_membranacea_Hawaii")

temp <- subset(accessionData, tiplabel %in% pepRAXML_rooted$tip.label)

sort(unique(temp$newtiplabel))

# Define colors
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
plot(1:8, col = c(col1, col2, col3, col4, col5, col6, col7, col8), cex = 3, pch = 16)

legendDF <- data.frame(cat=sort(unique(accessionData$HigherGeography)), x = 1:9, y = 1.9)
testLegend <- ggplot(data = legendDF) + geom_point(aes(y = y, x = x, color = cat), size = 3) + scale_colour_manual(values = colBiogeog) + theme(legend.key = element_blank(), legend.title = element_blank())

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

plotLegend <- g_legend(testLegend)
ggsave(plotLegend, filename = file.path(phylo.dir, "biogeogLegend.pdf"), device = cairo_pdf, width = 2, height = 4)

# For this method to work, the FIRST column must be matched with the tip labels
phyloSupport <- ggtree(pepRAXML_rooted, size = 0.2, ladderize = TRUE) %<+%
  accessionData[c("tiplabel","newtiplabel", "HigherGeography")] +
  #geom_text2(aes(label=label, subset=!is.na(as.numeric(label)) & as.numeric(label)>80 ), size = 3, hjust = -.2) + # plots bootstraps
  geom_nodepoint(aes(subset=!is.na(as.numeric(label)) & as.numeric(label)>80), color = "grey30", alpha = 0.5, size= 3) + # adds a dot for nodes > 80 bootstrap
  #geom_nodepoint(aes(subset=!is.na(as.numeric(label)) & as.numeric(label)>70 & as.numeric(label)<=80), color = "grey30", alpha = 0.5, size= 2) + 
  geom_hilight(node = 214, fill = "grey90", extendto = 0.035) +
  geom_hilight(node = 132, fill = "grey90", extendto = 0.035) +
  geom_tiplab(aes(label=newtiplabel, color = HigherGeography), size = 1.5) +
  #xlim(0, 0.045) +
  xlim(0, 0.045) +
  #geom_cladehi(node = 214, label = "Hawaiian endemics A") +
  #geom_cladehi(node = 132, label = "Hawaiian endemics B") +
  scale_color_manual(values = colBiogeog) +
  geom_treescale(x = 0, y = 0, offset = -2) +
  theme(panel.background = element_rect(fill = "transparent"))

phyloNode <- ggtree(pepRAXML_rooted, size = 0.2, ladderize = TRUE) + 
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 1)

#pepPhylogram <- arrangeGrob(phyloSupport, g_legend(testLegend), nrow = 1, ncol = 2, widths = c(8, 2)) 

ggsave(phyloSupport, width = 8, height = 12, filename = file.path(phylo.dir, "pepPhylo_support.pdf"), device = cairo_pdf)
#ggsave(pepPhylogram, width = 8, height = 12, filename = file.path(phylo.dir, "pepPhylogram.pdf"), device = cairo_pdf)

ggsave(phyloNode, width = 8, height = 12, filename = file.path(phylo.dir, "pepPhylo_node.pdf"), device = cairo_pdf)

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







# Solution: https://stackoverflow.com/questions/46918358/ggplot2-add-geologic-time-axis-annotation-filled-rects-with-names-between-axis
library(magrittr)
library(dplyr)
library(gtable)
GTS <- read.csv("https://raw.githubusercontent.com/japhir/stratPlot/master/GTS_colours.csv") %>%
  mutate(col = rgb(R, G, B, maxColorValue = 255)) %>%
  mutate(mean = (end - start) / 2 + start)

y.range <- c(0, 60)
filter.GTS <- GTS %>% 
  mutate(type = factor(type, levels = c("Eon", "Era", "Period", "Epoch", "Age"))) %>%
    # optional: plot only subset of bar types
  filter(type %in% c("Epoch")) %>%
  filter(end >= min(y.range),
         start <= max(y.range)) %>%
  rowwise() %>%
  mutate(start = max(start, min(y.range)),
         end = min(end, max(y.range))) %>%
  ungroup() %>%
  mutate(height = (end - start) / (max(end) - min(start))) %>%
  
  # optional: suppress names for short bars (e.g. Bartonian Age)
  # if bar height < 10% of plot height
  # mutate(name = ifelse(height <= 0.1, "", name)) %>%
  
  select(name, type, col, height) %>%
  arrange(type)

unique.types <- unique(filter.GTS$type) %>% as.character()

# create empty gtable
gt <- gtable(widths = rep(unit(1, "null"), 
                          times = length(unique.types)),
             heights = unit(1, "null"))

# fill gtable with individual table grobs for each type of geologic time

period.df <- filter(filter.GTS, type == unique.types)
tt <- tableGrob(d = select(period.df, name),
                heights = unit(period.df$height, "null"),
                widths = unit(1, "null"),
                cols= NULL, rows = NULL,
                theme = ttheme_minimal(core = list(bg_params = list(fill = "white",
                                               col = "black")))
                )

gt <- gtable_add_grob(x = gt,
                      grobs = tt,
                      t = 1, l = i)
grid.draw(gt)

for(i in 1:length(period.df$name)){
  tt <- tableGrob(d = period.df,
                  cols = NULL, rows = NULL,
                  heights = unit(period.df$height, "null"),
                  widths = unit(1, "null"),
                  theme = ttheme_minimal(
                    core = list(bg_params = list(fill = "white",#period.df$col,
                                                 col = "black"))#,
                  )            #fg_params = list(rot = 90))
  )
  
}

gt <- gtable_add_grob(x = gt, grobs = tt, t = 1, l = 1)
  
testPhylo_node <- ggtree(pepRAXML_rooted) + geom_tiplab(size = 1.2) + geom_text(aes(label=node))

ggsave(testPhylo_node, width = 20, height = 20, filename = "~/Desktop/pepNode.pdf")

## COUNT

readDepth <- read.csv(file.path(main.dir, "readDepthSummary.csv"))
accessionData <- merge(readDepth, accessionData, by.x = "sample_ID", by.y = "SampleID" )
names(accessionData)
temp <- subset(accessionData, tiplabel %in% pepRAXML_rooted$tip.label)

quantile(temp$meanCoverage)
sd(temp$meanCoverage)
range(temp$meanCoverage)

mean(temp$signif.propAmbiguous..3.)
table(temp$HigherGeography)
nrow(subset(temp, HigherGeography %in% c("Asia", "Africa", "S + C America")))

## FIT BIOGEOBEARS ==============
pepRAXML