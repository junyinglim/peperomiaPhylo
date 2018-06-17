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

clockAnalysisList = ["Piper_kadsura", "PEZ-164_bwa", "PEZ-204_bwa", "PEZ-247_bwa", 
					 "PEZ-183_bwa", "PEZ-174_bwa", "PEZ-179_bwa", "PEZ-172_bwa", "PEZ-269_bwa", "PEZ-213_bwa", "PEZ-254_bwa", "PEZ-180_bwa",
					 "PEZ-198_bwa", "PEZ-193_bwa", "PEZ-197_bwa", "PEZ-218_bwa", "PEZ-224_bwa", "JYL-40_bwa", "PEZ-287_bwa", "JYL-71_bwa", "PEZ-184_bwa",
					 "PEZ-251_bwa", "PEZ-304_bwa", "PEZ-215_bwa", "JYL-51_bwa", "PEZ-248_bwa", "PEZ-253_bwa",
					 "JYL-37_bwa", "JYL-55_bwa", "PEZ-219_bwa", "PEZ-221_bwa", "JYL-42_bwa", "JYL-53_bwa", "PEZ-202_bwa", "PEZ-186_bwa", "PEZ-185_bwa", "JYL-56_bwa", "PEZ-136_bwa", "JYL-68_bwa", "JYL-61_bwa", "PEZ-226_bwa", "PEZ-203_bwa",
					 "PEZ-296_bwa", "PEZ-293_bwa", "PEZ-288_bwa", "PEZ-105_bwa", "PEZ-209_bwa", "PEZ-302_bwa", "PEZ-297_bwa", "PEZ-102_bwa", "PEZ-210_bwa", "PEZ-169_bwa", "PEZ-285_bwa", "PEZ-108_bwa", "PEZ-250_bwa", "PEZ-246_bwa", "PEZ-177_bwa",
					 "PEZ-281_bwa", "PEZ-278_bwa", "PEZ-294_bwa", "PEZ-166_bwa", "PEZ-168_bwa"]
# 1 = Piper
# 2 = S.Amer + tetraphylla
# 3 = Hawaii_A + Marq, (198=eekana, 193=membr. 197=kok., 218=pallida, 224=ligust., 40=remyi, 287=oliveri, 71=cookiana, 184=globulanthera)  
# 4 = S. Amer + blandas (251=coroicoensis, 304=rubella, 215=blanda Marq, 51=blanda Haw, 248=blanda Afrc, 253=fernandopoiana Afrc)
# 5 = Hawaii_B (37=maui, 55=sandw, 219=lat, 221=hirt, 42=alter, 53=oahu, 202=expall, 186=obov, 185=kipah, 56=hesper, 136=hypoleuca, 68=macrae, 61=ellipti, 226=subpet, 203=rockii)
# 6 = West + S. Pacific (296=enervis, 293=hunteriana, 288=tooviana, 105=societatis, 209=ponapensis, 302=bonin, 297=bellend, 102=palauensis, 210=glassmanii, 169=reineckei, 285=grantii, 108=hombronii, 250=caledonica, 246=caldenoiasp., 177=urviellana
# 7 = Fiji + Samoa (281+278=fiji, 294=adamsonii, 168=samoensis)

# Total=63 taxa

alignmentClockAnalysisList = [record for record in concatenatedAlignment if record.id in clockAnalysisList]
concatenatedAlignmentClockAnalysis = MultipleSeqAlignment(alignmentClockAnalysisList)

AlignIO.write(concatenatedAlignmentClockAnalysis, os.path.join(concat_dir, "concatAlign_clock_160618.phy"), "phylip-relaxed")

