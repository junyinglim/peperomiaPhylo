#!/bin/bash

## Extract Annotations using the reference genome 

# Packages --------------------
from Bio import SeqIO, AlignIO
import os
import subprocess
import re
import pandas as pd

# Directories for local --------------------
mafft="/usr/local/bin/mafft"
blat="/Users/junyinglim/Dropbox/Projects/Programs/blat/blat/blat"
assembled_dir="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/data/chloroplast_bwa_assembled" # path with consensus assemblies

# Import annotated genome --------------------
annotatedGenome=SeqIO.read("/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/referenceGenome_KT223569.1.gb", "genbank")

# Generate a unique list of gene ids --------------------
proteinCodingGenes = [feature for feature in annotatedGenome.features if feature.type == "CDS"]
proteinCodingGenes_names = [cds.qualifiers["gene"][0] for cds in proteinCodingGenes]
proteinCodingGenes_names = list(set(proteinCodingGenes_names))

# if there is more than 1 feature of that name, pick the first one    
proteinCodingGenesFeatures = []
for i in proteinCodingGenes_names:
    cds = [feature for feature in proteinCodingGenes if feature.qualifiers["gene"][0] == i]
    proteinCodingGenesFeatures.append(cds[0])

# if there are more than 2 parts, stitch them together
proteinCodingGenesSeq = []
for j in proteinCodingGenesFeatures:
    # if it is a compound location
    if len(j.location.parts) > 1: 
        seq = annotatedGenome[0:0] # empty SeqRecord to append
        for part in j.location.parts:
            seq = seq + annotatedGenome[int(part.start):int(part.end)]
        proteinCodingGenesSeq.append(seq)
    else:
       proteinCodingGenesSeq.append(annotatedGenome[int(j.location.start):int(j.location.end)])

# rename fastas
for k in range(len(proteinCodingGenesSeq)):
    proteinCodingGenesSeq[k].id = proteinCodingGenes_names[k]
    proteinCodingGenesSeq[k].name = proteinCodingGenes_names[k]
    proteinCodingGenesSeq[k].description = proteinCodingGenes_names[k]


# Extract protein coding genes into a single fasta file
proteinCodingGenesSeqFile = os.path.join(assembled_dir, "proteinCodingGenes.fasta")

SeqIO.write(proteinCodingGenesSeq, proteinCodingGenesSeqFile, "fasta")

# Import assembled plastomes --------------------
assemblies = [x for x in os.listdir(assembled_dir) if "assembled.fasta" in x]
assembliesFiles = [os.path.join(assembled_dir, x) for x in assemblies]
assembliesSeq = [SeqIO.read(x, "fasta") for x in assembliesFiles]

for i in range(len(assemblies)):
    assembliesSeq[i].description =  re.sub("_assembled.fasta", "", assemblies[i])
    assembliesSeq[i].id =  re.sub("_assembled.fasta", "", assemblies[i])

assembliesSeq.append(annotatedGenome) # add the reference genome to the list for alignment

seqFile = os.path.join(assembled_dir, "pep_seq.fasta")
SeqIO.write(assembliesSeq, seqFile, "fasta")

# Align platomes using MAFFT --------------------
alignmentFile = os.path.join(assembled_dir, 'pep_seq_aligned.fasta')
res = subprocess.check_output( [mafft, '--auto', seqFile])

f = open(alignmentFile, 'wb+')
f.write(res)
f.close


# Extract reference sequence from alignment
pep_alignment = AlignIO.read(alignmentFile, "fasta")
alignedReference =[seq for seq in pep_alignment if seq.id == "KT223569.1"][0]

alignedReferenceFile = os.path.join(assembled_dir, "alignedReference.fasta")
SeqIO.write(alignedReference, alignedReferenceFile, "fasta")

# BLAT back to reference sequence --------------------
blatOutputFile = os.path.join(assembled_dir, 'blatOutput.fsl')
subprocess.Popen( [blat, alignedReferenceFile, proteinCodingGenesSeqFile, blatOutputFile, '-noHead', '-minIdentity=100'] )

# Extract out from alignments
# alignPositions = pd.read_table(blatOutputFile,
#     names = ["matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert",
#     "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", "qStart", "qEnd",
#     "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts"])

# alignPositions = alignPositions[alignPositions['strand'] == '+'] # only use those on one strand

# Maybe I should not bother with CDS, and just deal with the full gene sequence