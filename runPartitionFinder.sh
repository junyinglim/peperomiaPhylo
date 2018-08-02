#!/bin/bash

partitionfinder="/Users/junyinglim/Desktop/partitionfinder-2.1.1/PartitionFinder.py"
#dir="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/raxml_180718/partition_finder_GTRGAMMAINV"
dir="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/pepPhyloRuns/raxml_180718/partition_finder_GTRGAMMA"
python $partitionfinder $dir/bothConcat_180718_partition.phy $dir/partition_finder.cfg --raxml --quick --processes=2


#python $partitionfinder $dir/bothConcat_180718.phy $dir/partitionFinder_GTRIG_both.cfg --raxml --quick --processes=2


