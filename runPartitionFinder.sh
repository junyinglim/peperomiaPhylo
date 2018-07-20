#!/bin/bash

partitionfinder="/Users/junyinglim/Desktop/partitionfinder-2.1.1/PartitionFinder.py"
dir="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/raxml_180718/"
python $partitionfinder $dir/codingConcat_180718.phy $dir/partition_finder.cfg --raxml --quick --processes=2