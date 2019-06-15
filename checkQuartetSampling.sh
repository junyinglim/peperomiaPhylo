# Quantifying branch support values for a phylogenetic tree using a quartet sampling method
# Accessed from: https://github.com/FePhyFoFum/quartetsampling
# Reference: James B Pease, Joseph W Brown, Joseph F Walker, Cody E Hinchliff, Stephen A Smith. 2018. Quartet Sampling distinguishes lack of support from conflicting support in the green plant tree of life. American Journal of Botany. doi:10.1002/ajb2.1016

quartet="/Users/junyinglim/Dropbox/Projects/Programs/quartetsampling-master/pysrc"

seq_dir="/Users/junyinglim/Dropbox/Peperomia/Pacific/genomeskimming_data/pepPhyloRuns/2018-09-04"
alignment="$seq_dir/bothTrim_2018-09-04_QS.phy"
partfile="$seq_dir/both_trim_partition_2018-09-04.txt"
treefile="$seq_dir/bothTrim_2018-09-04_IQTREE_QS.nwk"


python3 $quartet/quartet_sampling.py -t $treefile -a $alignment -q $partfile -N 100 -T 4 -L 2 -d nuc --raxml-model GTRGAMMA -X raxmlHPC-PTHREADS-SSE3 -o $seq_dir/qs/ -e $seq_dir/qs/QuartetSampling
