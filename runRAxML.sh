#!/bin/bash
# Build index file

raxmlHPC-PTHREADS-SSE3 -s concatAlign_310818.phy -n pep -m GTRGAMMAI -x 12345 -p 12345 -N 1000 -f a -T 4
raxmlHPC-PTHREADS-SSE3 -s concatAlign_310818.phy -n pep -m GTRGAMMAI -x 12345 -p 12345 -N 1000 -f a -T 4 -q pepPhylo_partition.txt
