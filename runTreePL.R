## Run TreePL on posterior sample of 

fbdRes <- read.delim("~/Desktop/fossilBD/data/fbd_uexp_test8.log")
rootAge <- fbdRes$crown_age_piperaceae

main.dir <- "/Users/junyinglim/Desktop/testTreePL"

rootAgePost <- sample(rootAge, size = 100)

generateTreePLconfig <- function(x, iter){
  treefile = paste0("treefile = ", "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/2018-09-04/iqtree/output.treefile_rooted_treePL")
  outfile = "outfile = treePL.tre"
  numsites = "numsites = 114191"
  mrca = "mrca = PIPERACEAE JYL-03_bwa Piper_kadsura"
  min = paste("min = PIPERACEAE", x)
  max = paste("max = PIPERACEAE", x)
  cvoutfile = "cvoutfile = cv.out"
  cv = "cv"
  thorough = "thorough"
  cvstart = "cvstart = 0.1"
  cvstop = "cvstop = 1000"
  fileConn <- file(file.path(main.dir, "treePL.config"))
  writeLines(c(treefile, outfile, numsites, mrca, min, max, cvoutfile, cv, thorough, cvstart, cvstop), fileConn)
}

library(subprocess)

for(i in rootAgePost){
  i = rootAgePost[1]
  generateTreePLconfig(i)
  #system2("/usr/local/bin/treepl", file.path(main.dir, "treePL.config"), wait = TRUE)
  spawn_process("/usr/local/bin/treepl", file.path(main.dir, "treePL.config"), workdir = main.dir)
}

## FOSSIKLIZED BIRTH DEATH TREES
library(ape)
library(ggtree)

fbdPostMAPTree <- read.beast("~/Desktop/fossilBD/data/fbd_uexp_test8_map.tree")
fossilTaxa <- c("Hexagyne_philippiana", "Lactoripollenites_africanus","Piper_margaritae", "Piper_bartlingianum", "Saururus_aquilae", "Saururus_tuckerae", "Saururus_bliobatus", "Saururus_stoobensis", "Houttuynia_bavarica", "Saururopsis_niponensis")


fbdPostMCCTreeExtant <- drop.tip(fbdPostMCCTree, tip = fossilTaxa)
