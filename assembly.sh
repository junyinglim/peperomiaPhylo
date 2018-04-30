#!/bin/bash

# Dependencies:
#   samtools, bowtie2, trimmomatic, gunzip
#   bcftools (tools for variant calling)

# Directories for the server
# DIR="/home/junyinglim/peperomiaFiles"
# trimmomatic="/home/junyinglim/Trimmomatic-0.36"
# indexedGenome="/home/junyinglim/referenceGenome/piperchloroplast"

# Directories for local
DIR="/Users/junyinglim/Desktop/testtest"
trimmomatic="/Users/junyinglim/Dropbox/Projects/Programs/Trimmomatic-0.36"
indexedGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/piperchloroplast"
referenceGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/referenceGenome_KT223569.1.fasta"

for s in $DIR/*_L008_R1*.fastq.gz
do
    s=${s##*/}
    index=${s%_S*}

    # echo "Uncompressing index: $index"
    # gunzip $DIR/${index}_*R1*.fastq.gz -k -d
    # gunzip $DIR/${index}_*R2*.fastq.gz -k -d
    
    # echo "Filtering adapter sequences: $index"
    # java -jar $trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 $DIR/${index}_*R1*.fastq $DIR/${index}_*R2*.fastq $DIR/${index}_trimmed_P1.fastq $DIR/${index}_trimmed_U1.fastq $DIR/${index}_trimmed_P2.fastq $DIR/${index}_trimmed_U2.fastq ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE.fa:1:30:10 SLIDINGWINDOW:10:20 MINLEN:40

    # echo "Map to chloroplast reference using bowtie2: $index"
    # bowtie2 -1 $DIR/${index}_trimmed_P1.fastq -2 $DIR/${index}_trimmed_P2.fastq -x $indexedGenome -S $DIR/${index}_mapped.sam

    # echo "Convert sam to bam: $index" # Converts from a human-readable format to a binary format
    # samtools view -bS $DIR/${index}_mapped.sam > $DIR/${index}_mapped.bam

    # echo "Sorting reads: $index" # Sorts reads by their aligned position to reference
    # samtools sort $DIR/${index}_mapped.bam > $DIR/${index}_sorted.bam

    # echo "Determining read depth: $index" # Outputs sequencing depth in a tab delimited format
    # samtools depth $DIR/${index}_sorted.bam > $DIR/${index}_read_depth.tsv

    # echo "Filtering by phred quality 20: $index" # Creates a pileup file; summarizes the base calls of aligned reads
    # samtools mpileup -Q 20 -Agf $referenceGenome $DIR/${index}_sorted.bam > $DIR/${index}.mpilup

    # echo "Generating consensus genotypes..."
    # # 
    # bcftools call --ploidy 1 -c $DIR/${index}.mpilup > $DIR/${index}.vcf  #Assuming ploidy is 1 for chloroplast 

    echo "Filtering for read depth >= 10..."
    vcftools --vcf $DIR/${index}.vcf --out $DIR/${index}_filter.vcf --minDP 10 --recode --recode-INFO-all

    echo "Convert vcf file back into a fasta"
    vcfutils.pl vcf2fq ${index}.filter.vcf > ${index}_assembled.fastq
    seqtk seq -A ${index}_assembled.fastq > ${index}_assembled.fasta


    #echo "Cleaning up for now.."

        #mkdir ${output_dir}${index}
        #mv ${index}.fasta ${index}_read_depth.tsv ${output_dir}${index}
        

    bedtools # convert from bam to fasta
    # find ORF

    # will the consensus genome have gaps?

done