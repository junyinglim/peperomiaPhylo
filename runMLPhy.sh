#!/bin/bash
# Build index file



# -n output file name
# -m substitution model
# -x turns on fast bootstrapping, and sets seed
# -N = number of rapid bootstrap trees
# -d start Raxml search with a random tree (if not specified, the default is a maximum parsimony tree)


#raxmlHPC-PTHREADS-SSE3 -s concatAlign_160818.phy -n 160818 -m GTRGAMMA -x 12345 -p 12345 -N 1000 -f a -T 4

DIR="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/raxml_180718"
iqtree -s $DIR/bothConcat_180718.phy -alrt 1000 -bb 1000 -spp $DIR/both_partition_180718.txt -m MFP+MERGE+R -nt 4

# bb = number of ultrafast boostraps
# alrt = nrach tests
# spp = edge-linked proportional partition model
# -m MFP+MERGE = find best partition scheme (including the free rate heterogeneity model), followed by tree inference
# -nt = 4 CPU cores

