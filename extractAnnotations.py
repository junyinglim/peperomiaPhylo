#!/bin/bash

## Extract Annotations using the reference genome
from Bio import SeqIO, AlignIO
import os
import subprocess
import re

# Directories for local
DIR="/Users/junyinglim/Desktop/testtest"
mafft="/usr/local/bin/mafft"

# Import annotated genome
annotatedGenome=SeqIO.read("/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/referenceGenome_KT223569.1.gb", "genbank")
proteinCodingGenes = [feature for feature in annotatedGenome.features if feature.type == "CDS"]

proteinCodingGenes_names = [cds.qualifiers["gene"] for cds in proteinCodingGenes]
proteinCodingGene_start = [int(cds.location.start) for cds in proteinCodingGenes]
proteinCodingGene_end = [int(cds.location.end) for cds in proteinCodingGenes]

# Extract protein coding genes into a single fasta file



# Import assembled plastomes
assembledDIR="/Users/junyinglim/Desktop/assembled"
assemblies = [x for x in os.listdir(assembledDIR) if "assembled.fasta" in x]
assembliesFiles = [os.path.join(assembledDIR, x) for x in assemblies]

assembliesSeq = [SeqIO.read(x, "fasta") for x in assembliesFiles]

for i in range(len(assemblies)):
    assembliesSeq[i].description =  re.sub("_assembled.fasta", "", assemblies[i])
    assembliesSeq[i].id =  re.sub("_assembled.fasta", "", assemblies[i])

assembliesSeq.append(annotatedGenome) # add the reference genome into the list o

seqfile = os.path.join(assembledDIR, "pep_seq.fasta")
SeqIO.write(assembliesSeq, seqfile, "fasta")

# Align platomes using MAFFT
outfile = os.path.join(assembledDIR, 'pep_seq_aligned.fasta')
res = subprocess.check_output( [mafft, '--auto', seqfile])
with file(outfile, 'w') as f:
    f.write(res)

# Extract reference sequence
pep_alignment = AlignIO.read(outfile, "fasta")
alignedReference =[seq for seq in pep_alignment if seq.id == "KT223569.1"][0]

# BLAT back to reference sequence


# Extract reference sequence (with gaps)
# Map genes onto the alignment
# Get the locations on reference sequence
# Use Stefan's code (or write custom code) to extract those genes individually
