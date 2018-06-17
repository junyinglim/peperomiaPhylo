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

alignmentList = [record for record in concatenatedAlignment if record.id not in ["PEZ-217_bwa",
		"PEZ-231_bwa", "PEZ-232_bwa", "PEZ-233_bwa", "PEZ-234_bwa", "PEZ-235_bwa", "PEZ-236_bwa",
		"PEZ-237_bwa", "PEZ-238_bwa", "PEZ-239_bwa", "PEZ-241_bwa", "PEZ-242_bwa", "PEZ-243_bwa",
		"PEZ-244_bwa", "PEZ-245_bwa"]]
# remove Piper (204 + 247 + 164) ("Piper_kadsura", "PEZ-164_bwa", "PEZ-204_bwa", "PEZ-247_bwa")
# 234, 235, 237, 241 don't have locality information (P. blandas and P. polystachya)
# 217 (weird Pep blanda) may be an error
# 231, 232, 233 (austral blanda, no voucher?)
# 236 (pitcairns, no voucher?)
# 239, 241, 242, 243, 244 (Juan F.); PEZ-241; PEZ-242; PEZ-243 (Juan F.), PEZ-244 (Juan F.) don't have a living collections number

concatenatedAlignment_subset = MultipleSeqAlignment(alignmentList)

AlignIO.write(concatenatedAlignment_subset, os.path.join(concat_dir, "concatAlign_160618.phy"), "phylip-relaxed")


#alignmentList = [record for record in concatenatedAlignment if record.id in ["PEZ-179_bwa",]]
