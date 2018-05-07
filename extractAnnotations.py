#!/bin/bash

## Extract Annotations using the reference genome
from Bio import SeqIO

# Directories for local
DIR="/Users/junyinglim/Desktop/testtest"
indexedGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/piperchloroplast"
referenceGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/referenceGenome_KT223569.1.fasta"
mafft="/usr/local/bin/mafft"


# Import annotated genome
annotatedGenome=SeqIO.read("/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/referenceGenome_KT223569.1.gb", "genbank")
proteinCodingGenes = [feature for feature in annotatedGenome.features if feature.type == "CDS"]

proteinCodingGenes_names = [cds.qualifiers["gene"] for cds in proteinCodingGenes]
proteinCodingGene_start = [int(cds.location.start) for cds in proteinCodingGenes]
proteinCodingGene_end = [int(cds.location.end) for cds in proteinCodingGenes]

# Sequence alignment of assembled plastomes




# Next extract the genes out into individual fastas
# MAFFT