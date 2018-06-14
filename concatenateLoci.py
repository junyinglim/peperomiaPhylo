#!/bin/bash

## Extract Annotations using the reference genome 

# Packages --------------------
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os
from seqTools import genPartition, concatenate # custom functions

# Directories for local --------------------
alignments_dir ="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/bwa_alignments"
concat_dir = "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data"

# Import locus alignments --------------------
filenames = os.listdir(alignments_dir)
fastas=[file for file in filenames if file.endswith(".fasta")]

fastas=[fasta for fasta in fastas if fasta != "ycf1_align.fasta"] # remove ycf1

alignments=[AlignIO.read(os.path.join(alignments_dir, fasta), "fasta") for fasta in fastas]
genPartition(alignments, filename = os.path.join(concat_dir, "pepPhylo_partition.txt"))

# Concatenate locus alignments -------------------- 
concatenatedAlignment = concatenate(alignments)

alignmentList = [record for record in concatenatedAlignment if record.id not in ["Piper_kadsura", "PEZ-250_bwa", "PEZ-251_bwa", "PEZ-217_bwa", "PEZ-164_bwa", "PEZ-204_bwa", "PEZ-247_bwa"]] 
# remove Piper (204 + 247 + 164), 250+251 (low coverage), 217 (mis-label)

concatenatedAlignment_subset = MultipleSeqAlignment(alignmentList)

AlignIO.write(concatenatedAlignment_subset, os.path.join(concat_dir, "concatAlign_130618.phy"), "phylip-relaxed")

# Generate partition parameters -------------------- 
