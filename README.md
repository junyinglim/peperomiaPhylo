# peperomiaPhylo


## Pipeline
* Adapter trimming using `Trimmomatic`
* Index the refence genome using `bowtie2` (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome)
* Mapping to reference genome using `bowtie2`


* Evaluating coverage of chloroplast genomes
* Extract out chloroplast genes using annotations from reference
* Run a phylogeny using chloroplast genes


## Ideas
* Streamline code (a few unnecessary steps?)
    
* Find ORF vs. using annotation codes to extract taxa
* Will the consensus genome have gaps?