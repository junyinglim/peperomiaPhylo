#!/bin/bash
# Build index file



# raxmlHPC-PTHREADS-SSE3 -s concatAlign_310818.phy -n pep -m GTRGAMMAI -x 12345 -p 12345 -N 1000 -f a -T 4
# raxmlHPC-PTHREADS-SSE3 -s concatAlign_310818.phy -n pep -m GTRGAMMAI -x 12345 -p 12345 -N 1000 -f a -T 4 -q pepPhylo_partition.txt


# Ribosome run
#raxmlHPC-PTHREADS-SSE3 -s ribosome_bwa_align.phy -n pepRibo -m GTRGAMMA -x 12345 -p 12345 -N 1000 -f a -T 4
#raxmlHPC-PTHREADS-SSE3 -s ribosome_bwa_align.phy -n pepRiboInv -m GTRGAMMAI -x 12345 -p 12345 -N 1000 -f a -T 4

# -n output file name
# -m substitution model
# -x turns on fast bootstrapping, and sets seed
# -N = number of rapid bootstrap trees
# -d start Raxml search with a random tree (if not specified, the default is a maximum parsimony tree)


raxmlHPC-PTHREADS-SSE3 -s concatAlign_160818.phy -n 160818 -m GTRGAMMA -x 12345 -p 12345 -N 1000 -f a -T 4