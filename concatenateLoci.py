#!/bin/bash

## Extract Annotations using the reference genome 

# Packages --------------------
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os
from seqTools import genPartition, concatenate # custom functions

from Bio.AlignIO import PhylipIO
import subprocess


import datetime

now = datetime.datetime.now()
stem = now.strftime("%Y-%m-%d")

# Directories for local --------------------
align_dir = "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/bwa_alignments/"
coding_dir = os.path.join(align_dir, "coding")
noncoding_dir =os.path.join(align_dir, "noncoding")

coding_trim_dir = os.path.join(align_dir, "coding_trim")
noncoding_trim_dir = os.path.join(align_dir, "noncoding_trim")

concat_dir = "/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data"

sed = "/usr/bin/sed"
trimal = "/Users/junyinglim/Dropbox/Projects/Programs/trimAl-v1.4/source/trimal"



# Import locus alignments --------------------
codingfastas=[file for file in os.listdir(coding_dir) if file.endswith(".fasta")]
noncodingfastas=[file for file in os.listdir(noncoding_dir) if file.endswith(".fasta")]

codingfastas=[fasta for fasta in codingfastas if fasta != "ycf1_align.fasta"] # remove ycf1

noncodingexclude = ["atpA-atpF_align.fasta", "atpF-atpH_align.fasta", "ccsA-trnL_align.fasta",
					"ndhB-rps7_align.fasta", "ndhD-ccsA_align.fasta", "ndhG-ndhE_align.fasta",
					"petD-rpoA_align.fasta", "psbE-petL_align.fasta", "rpl32-ndhF_align.fasta",
					"rps3-rpl22_align.fasta", "rps16-trnQ-psbK_align.fasta", "rps19-rpl2_align.fasta",
					"trnE-trnT_align.fasta", "trnN-ycf1_align.fasta", "ycf1-rps15_align.fasta"]
noncodingfastas=[fasta for fasta in noncodingfastas if fasta not in noncodingexclude]


codingAlignments=[AlignIO.read(os.path.join(coding_dir, fasta), "fasta") for fasta in codingfastas]
noncodingAlignments=[AlignIO.read(os.path.join(noncoding_dir, fasta), "fasta") for fasta in noncodingfastas]
bothAlignments = codingAlignments + noncodingAlignments

genPartition(codingAlignments,
	filename = os.path.join(concat_dir, "coding_partition_" + stem + ".txt"))
genPartition(noncodingAlignments,
	filename = os.path.join(concat_dir, "noncoding_partition_" + stem + ".txt"))
genPartition(bothAlignments,
	filename = os.path.join(concat_dir, "both_partition_" + stem + ".txt"))


# Concatenate locus alignments -------------------- 
concatenatedAlignment_coding = concatenate(codingAlignments)
concatenatedAlignment_noncoding = concatenate(noncodingAlignments)
concatenatedAlignment_both = concatenate(bothAlignments)

excludedSpecimens = ["PEZ-231_bwa", "PEZ-233_bwa", "PEZ-234_bwa", "PEZ-235_bwa", "PEZ-237_bwa", "PEZ-244_bwa", "PEZ-245_bwa"]
# 231,233,234,244, 245 = no more living collection
# 235 = voucher but no locality
# 237 = no locality

concatenatedAlignment_coding_subset = [record for record in concatenatedAlignment_coding if record.id not in excludedSpecimens]
concatenatedAlignment_noncoding_subset = [record for record in concatenatedAlignment_noncoding if record.id not in excludedSpecimens]
concatenatedAlignment_both_subset = [record for record in concatenatedAlignment_both if record.id not in excludedSpecimens]

concatenatedAlignment_coding_subset = MultipleSeqAlignment(concatenatedAlignment_coding_subset)

concatenatedAlignment_noncoding_subset = MultipleSeqAlignment(concatenatedAlignment_noncoding_subset)

concatenatedAlignment_both_subset = MultipleSeqAlignment(concatenatedAlignment_both_subset)

AlignIO.write(concatenatedAlignment_coding_subset,
	os.path.join(concat_dir, "codingConcat_"+stem+".phy"), "phylip-relaxed")
AlignIO.write(concatenatedAlignment_noncoding_subset,
	os.path.join(concat_dir, "noncodingConcat_"+stem+".phy"), "phylip-relaxed")
AlignIO.write(concatenatedAlignment_both_subset,
	os.path.join(concat_dir, "bothConcat_"+stem+".phy"), "phylip-relaxed")


# Triming using Trimal --------------------
[subprocess.Popen( [sed, "s/n/-/g", os.path.join(coding_dir, fasta)], stdout=open(os.path.join(coding_trim_dir, fasta), "w") ) for fasta in codingfastas]

[subprocess.Popen( [sed, "s/n/-/g", os.path.join(noncoding_dir, fasta)], stdout=open(os.path.join(noncoding_trim_dir, fasta), "w") ) for fasta in noncodingfastas]


[subprocess.Popen( [trimal, "-in", os.path.join(coding_trim_dir, fasta), "-out", os.path.join(coding_trim_dir, fasta), "-gappyout", "-fasta"] ) for fasta in codingfastas]
[subprocess.Popen( [trimal, "-in", os.path.join(noncoding_trim_dir, fasta), "-out", os.path.join(noncoding_trim_dir, fasta), "-gappyout", "-fasta"] ) for fasta in noncodingfastas]

codingTrimAlignments=[AlignIO.read(os.path.join(coding_trim_dir, fasta), "fasta") for fasta in codingfastas]
noncodingTrimAlignments=[AlignIO.read(os.path.join(noncoding_trim_dir, fasta), "fasta") for fasta in noncodingfastas]
bothTrimAlignments = codingTrimAlignments + noncodingTrimAlignments

genPartition(codingTrimAlignments,
	filename = os.path.join(concat_dir, "coding_trim_partition_" + stem + ".txt"))
genPartition(noncodingTrimAlignments,
	filename = os.path.join(concat_dir, "noncoding_trim_partition_" + stem + ".txt"))
genPartition(bothTrimAlignments,
	filename = os.path.join(concat_dir, "both_trim_partition_" + stem + ".txt"))



# Concatenate trimmed alignments --------------------
concatenatedTrimAlignment_coding = concatenate(codingTrimAlignments)
concatenatedTrimAlignment_noncoding = concatenate(noncodingTrimAlignments)
concatenatedTrimAlignment_both = concatenate(bothTrimAlignments)

removeSpecimens = ["PEZ-231_bwa", "PEZ-233_bwa", "PEZ-234_bwa", "PEZ-235_bwa", "PEZ-237_bwa", "PEZ-244_bwa", "PEZ-245_bwa"]

lowCoverageSpecimens = ["PEZ-251_bwa", "PEZ-250_bwa", "PEZ-209_bwa",
						"JYL-59_bwa", "PEZ-194_bwa", "PEZ-193_bwa",
						"PEZ-168_bwa"]
# <15x mean and median coverage
# 194 = kipahuluensis (no redund. for Kauai population)
# 193

# based on the ambiguity and gap trimmed alignment (and excluding the Pipers)
nonhomogeneousSpecimens = ["JYL-59_bwa", "PEZ-202_bwa",
						   "PEZ-205_bwa", "PEZ-209_bwa",
						   "PEZ-210_bwa"]

# JYL-59, (PEZ-164), PEZ-202, (204), 205, 209, 210, (247), (Piper_kadsura) [piper taxa in parentheses]
#205 = hirtipetiola (some redund.), 202 = expallescens (no redund.), 59 = macraeana (multiple redundancies)

excludedSpecimens = list(set().union(removeSpecimens,lowCoverageSpecimens, nonhomogeneousSpecimens))



concatenatedTrimAlignment_coding_subset = [record for record in concatenatedTrimAlignment_coding if record.id not in excludedSpecimens]
concatenatedTrimAlignment_noncoding_subset = [record for record in concatenatedTrimAlignment_noncoding if record.id not in excludedSpecimens]
concatenatedTrimAlignment_both_subset = [record for record in concatenatedTrimAlignment_both if record.id not in excludedSpecimens]

concatenatedTrimAlignment_coding_subset = MultipleSeqAlignment(concatenatedTrimAlignment_coding_subset)

concatenatedTrimAlignment_noncoding = MultipleSeqAlignment(concatenatedTrimAlignment_noncoding)

concatenatedTrimAlignment_both_subset = MultipleSeqAlignment(concatenatedTrimAlignment_both_subset)

AlignIO.write(concatenatedTrimAlignment_coding_subset,
	os.path.join(concat_dir, "codingTrim_"+stem+".phy"), "phylip-relaxed")

AlignIO.write(concatenatedTrimAlignment_noncoding,
	os.path.join(concat_dir, "noncodingTrim_"+stem+".phy"), "phylip-relaxed")

AlignIO.write(concatenatedTrimAlignment_both_subset,
	os.path.join(concat_dir, "bothTrim_"+stem+".phy"), "phylip-relaxed")




#alignmentList = [record for record in concatenatedAlignment if record.id in ["PEZ-179_bwa",]]

# clockAnalysisList = ["Piper_kadsura", "PEZ-164_bwa", "PEZ-204_bwa", "PEZ-247_bwa", 
# 					 "PEZ-183_bwa", "PEZ-174_bwa", "PEZ-179_bwa", "PEZ-172_bwa", "PEZ-269_bwa", "PEZ-213_bwa", "PEZ-254_bwa", "PEZ-180_bwa",
# 					 "PEZ-198_bwa", "PEZ-193_bwa", "PEZ-197_bwa", "PEZ-218_bwa", "PEZ-224_bwa", "JYL-40_bwa", "PEZ-287_bwa", "JYL-71_bwa", "PEZ-184_bwa",
# 					 "PEZ-251_bwa", "PEZ-304_bwa", "PEZ-215_bwa", "JYL-51_bwa", "PEZ-248_bwa", "PEZ-253_bwa",
# 					 "JYL-37_bwa", "JYL-55_bwa", "PEZ-219_bwa", "PEZ-221_bwa", "JYL-42_bwa", "JYL-53_bwa", "PEZ-202_bwa", "PEZ-186_bwa", "PEZ-185_bwa", "JYL-56_bwa", "PEZ-136_bwa", "JYL-68_bwa", "JYL-61_bwa", "PEZ-226_bwa", "PEZ-203_bwa",
# 					 "PEZ-296_bwa", "PEZ-293_bwa", "PEZ-288_bwa", "PEZ-105_bwa", "PEZ-209_bwa", "PEZ-302_bwa", "PEZ-297_bwa", "PEZ-102_bwa", "PEZ-210_bwa", "PEZ-169_bwa", "PEZ-285_bwa", "PEZ-108_bwa", "PEZ-250_bwa", "PEZ-246_bwa", "PEZ-177_bwa",
# 					 "PEZ-281_bwa", "PEZ-278_bwa", "PEZ-294_bwa", "PEZ-166_bwa", "PEZ-168_bwa"]
# # 1 = Piper
# # 2 = S.Amer + tetraphylla
# # 3 = Hawaii_A + Marq, (198=eekana, 193=membr. 197=kok., 218=pallida, 224=ligust., 40=remyi, 287=oliveri, 71=cookiana, 184=globulanthera)  
# # 4 = S. Amer + blandas (251=coroicoensis, 304=rubella, 215=blanda Marq, 51=blanda Haw, 248=blanda Afrc, 253=fernandopoiana Afrc)
# # 5 = Hawaii_B (37=maui, 55=sandw, 219=lat, 221=hirt, 42=alter, 53=oahu, 202=expall, 186=obov, 185=kipah, 56=hesper, 136=hypoleuca, 68=macrae, 61=ellipti, 226=subpet, 203=rockii)
# # 6 = West + S. Pacific (296=enervis, 293=hunteriana, 288=tooviana, 105=societatis, 209=ponapensis, 302=bonin, 297=bellend, 102=palauensis, 210=glassmanii, 169=reineckei, 285=grantii, 108=hombronii, 250=caledonica, 246=caldenoiasp., 177=urviellana
# # 7 = Fiji + Samoa (281+278=fiji, 294=adamsonii, 168=samoensis)

# # Total=63 taxa

# alignmentClockAnalysisList = [record for record in concatenatedAlignment if record.id in clockAnalysisList]
# concatenatedAlignmentClockAnalysis = MultipleSeqAlignment(alignmentClockAnalysisList)

# AlignIO.write(concatenatedAlignmentClockAnalysis, os.path.join(concat_dir, "concatAlign_clock_160618.phy"), "phylip-relaxed")

