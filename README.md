# peperomiaPhylo


## Pipeline
* Quality filtering of fastq files using fastq_quality_filter
* Index the refence genome using `bowtie2` (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome)
* Mapping to reference genome using `bowtie2`
* Evaluating coverage of chloroplast genomes
* Extract out chloroplast genes using annotations from reference
* Run a phylogeny using chloroplast genes
