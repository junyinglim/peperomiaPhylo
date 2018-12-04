import os
import pandas as pd
import re
import random
import subprocess

main_dir = "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/"
output_dir = os.path.join(main_dir, "output")

fbd_log_fnames = [x for x in os.listdir(output_dir) if ".log" in x]
fbd_log_fnames = [os.path.join(output_dir, x) for x in fbd_log_fnames]

fbd_logs = [pd.read_table(x) for x in fbd_log_fnames]

fbd_logs_noburnin = [x.loc[int(0.25*len(x.index)):len(x.index),] for x in fbd_logs ] # take the first 25% of iterations as burn-in
fbd_logs_noburnin = pd.concat(fbd_logs_noburnin) # rbind dataframes

meanDivTime = fbd_logs_noburnin["crownAge_pippep"].mean() # take the mean divergence time from the posterirop distribution

randomDivTime = random.sample(list(fbd_logs_noburnin["crownAge_pippep"]), 100) # randomly sample 100 from posterior distribution

#crownAgePipPep <- sample(fbd_logs_noburnin$crownAge_pippep, size = 100)

def generateTreePLconfig(x, iter):
	treefile = "treefile = " + "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/2018-09-04/iqtree/output.treefile_rooted_treePL"
	outfile = "outfile = treePL_" + str(iter) + ".out"
	numsites = "numsites = 114191"
	mrca = "mrca = PIPERACEAE JYL-03_bwa Piper_kadsura"
	min = "min = PIPERACEAE " + str(x)
	max = "max = PIPERACEAE " + str(x)
	cvoutfile = "cvoutfile = cv_" + str(iter) + ".out"
	cv = "cv"
	thorough = "thorough"
	cvstart = "cvstart = 0.1"
	cvstop = "cvstop = 1000"
	lines = [treefile, outfile, numsites, mrca, min, max, cvoutfile, cv, thorough, cvstart, cvstop]
	with open(os.path.join(output_dir, "treePL.config"), 'w') as f:
		f.write('\n'.join(lines))

generateTreePLconfig(x = meanDivTime, iter = "mean")
subprocess.Popen(["/usr/local/bin/treepl", os.path.join(output_dir, "treePL.config")])

# for i in range(len(randomDivTime)): 
# 	iter = i+1
# 	generateTreePLconfig(x = randomDivTime[i], iter = iter)
# 	subprocess.Popen(["/usr/local/bin/treepl", os.path.join(output_dir, "treePL.config")])


## FOSSIKLIZED BIRTH DEATH TREES
# library(ape)
# library(ggtree)
# 
# fbdPostMAPTree <- read.beast("~/Desktop/fossilBD/data/fbd_uexp_test8_map.tree")
# fossilTaxa <- c("Hexagyne_philippiana", "Lactoripollenites_africanus","Piper_margaritae", "Piper_bartlingianum", "Saururus_aquilae", "Saururus_tuckerae", "Saururus_bliobatus", "Saururus_stoobensis", "Houttuynia_bavarica", "Saururopsis_niponensis")
# 
# 
# fbdPostMCCTreeExtant <- drop.tip(fbdPostMCCTree, tip = fossilTaxa)
