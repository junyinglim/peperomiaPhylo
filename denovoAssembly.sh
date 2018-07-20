#!/bin/bash

# Dependencies:
#   samtools, bowtie2, trimmomatic, gunzip
#   bcftools (tools for variant calling)

# Directories for the server
# DIR="/home/junyinglim/peperomiaFiles"
# trimmomatic="/home/junyinglim/Trimmomatic-0.36"
# indexedGenome="/home/junyinglim/referenceGenome/piperchloroplast"
# referenceGenome="/home/junyinglim/referenceGenome/referenceGenome_KT223569.1.fasta"

# # Directories for local
DIR="/Users/junyinglim/Desktop/testtest"
trimmomatic="/Users/junyinglim/Dropbox/Projects/Programs/Trimmomatic-0.36"
indexedGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/piperchloroplast"
referenceGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/referenceGenome_KT223569.1.fasta"
#getorganelle="/Users/junyinglim/Dropbox/Projects/Programs/GetOrganelle/get_organelle_reads.py"
#spades="/Users/junyinglim/Dropbox/Projects/Programs/SPAdes-3.12.0-Darwin/bin"
fastplast="/Users/junyinglim/Desktop/Fast-Plast/fast-plast.pl"

for s in $DIR/*_L008_R1*.fastq.gz
do
    s=${s##*/}
    index=${s%_S*}

    # printf "\nUncompressing index: $index"
    # gunzip $DIR/${index}_*R1*.fastq.gz -k -d
    # gunzip $DIR/${index}_*R2*.fastq.gz -k -d
    
    # printf "\nFiltering adapter sequences: $index\n"
    # java -jar $trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 $DIR/${index}_*R1*.fastq $DIR/${index}_*R2*.fastq $DIR/${index}_trimmed_P1.fastq $DIR/${index}_trimmed_U1.fastq $DIR/${index}_trimmed_P2.fastq $DIR/${index}_trimmed_U2.fastq ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE.fa:1:30:10 SLIDINGWINDOW:10:20 MINLEN:40

    # printf "\nUsing GetOrganelle script to filter out organelle sequences\n"
    # python $getorganelle -1 $DIR/${index}_trimmed_P1.fastq -2 $DIR/${index}_trimmed_P2.fastq -s $referenceGenome -w 103 -o $DIR/${index}_chloroplast -R 5 -k 75,85,95,105 -P 300000 -F cp

    #printf "\n Testing spades"
    #python $spades/spades.py --pe1-1 $DIR/${index}_chloroplast/filtered_1_paired.fq --pe1-2 $DIR/${index}_chloroplast/filtered_2_paired.fq -t 2 -k 21,33,55,77,99,127 --careful -o $DIR/${index}_chloroplast/spades_output
    #get_organelle_reads.py -1 sample_1.fq -2 sample_2.fq -s cp_reference.fasta -w 103 -o chloroplast_output -R 5 -k 75,85,95,105 -P 1000000 -a mitochondria.fasta -J 3 -M 5

    #printf "\nUsing velvet to construct contigs: $index"
    #velveth $DIR/${index}_chloroplast/velvet 31 -shortPaired -fastq -separate $DIR/${index}_chloroplast/filtered_1_paired.fq $DIR/${index}_chloroplast/filtered_2_paired.fq
    #velvetg $DIR/${index}_chloroplast/velvet -amos_file yes -read_trkg yes

    printf "\nUsing Fast-plast to assemble genomes: $index\n"
    cd $DIR
    perl $fastplast -1 $DIR/${index}_*R1*.fastq -2 $DIR/${index}_*R2*.fastq -n ${index}_denovo --adapters Truseq --threads 4 --coverage_analysis --user_bowtie $indexedGenome 
    
done
#-s [$DIR/${index}_trimmed_U1.fastq,$DIR/${index}_trimmed_U2.fastq]


# Paths
#"/Users/junyinglim/Dropbox/Projects/Programs/Trimmomatic-0.36"
#"/Users/junyinglim/Dropbox/Projects/Programs/SPAdes-3.12.0-Darwin/"
#"/Users/junyinglim/Dropbox/Projects/Programs/sspace_basic"

# Had to go into the afin folder and make