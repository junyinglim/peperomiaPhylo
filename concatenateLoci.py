#!/bin/bash

## Extract Annotations using the reference genome 

# Packages --------------------
from Bio import AlignIO
import os
from seqTools import genPartition, concatenate # custom functions

# Directories for local --------------------
alignments_dir ="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/bwa_alignments"
concat_dir = "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data"

# Import locus alignments --------------------
filenames = os.listdir(alignments_dir)
fastas=[file for file in filenames if file.endswith(".fasta")]
alignments=[AlignIO.read(os.path.join(alignments_dir, fasta), "fasta") for fasta in fastas]

# Concatenate locus alignments -------------------- 
concatenatedAlignment = concatenate(alignments)
AlignIO.write(concatenatedAlignment, os.path.join(concat_dir, "concatAlign_310818.phy"), "phylip-relaxed")

# Generate partition parameters -------------------- 
genPartition(alignments, filename = os.path.join(concat_dir, "pepPhylo_partition.txt"))