#!/bin/bash

# Dependencies:
#   samtools, bowtie2, trimmomatic, gunzip
#   bcftools (tools for variant calling)

# Directories for the server
DIR="/home/junyinglim/peperomiaFiles"
trimmomatic="/home/junyinglim/Trimmomatic-0.36"
indexedGenome="/home/junyinglim/referenceGenome/piperchloroplast"
referenceGenome="/home/junyinglim/referenceGenome/referenceGenome_KT223569.1.fasta"

# # Directories for local
# DIR="/Users/junyinglim/Desktop/testtest"
# trimmomatic="/Users/junyinglim/Dropbox/Projects/Programs/Trimmomatic-0.36"
# indexedGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/piperchloroplast"
# referenceGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/referenceGenome_KT223569.1.fasta"

for s in $DIR/*_L008_R1*.fastq.gz
do
    s=${s##*/}
    index=${s%_S*}

    # printf "\nUncompressing index: $index"
    # gunzip $DIR/${index}_*R1*.fastq.gz -k -d
    # gunzip $DIR/${index}_*R2*.fastq.gz -k -d
    
    printf "\nFiltering adapter sequences: $index"
    java -jar $trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 $DIR/${index}_*R1*.fastq $DIR/${index}_*R2*.fastq $DIR/${index}_trimmed_P1.fastq $DIR/${index}_trimmed_U1.fastq $DIR/${index}_trimmed_P2.fastq $DIR/${index}_trimmed_U2.fastq ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE.fa:1:30:10 SLIDINGWINDOW:10:20 MINLEN:40

    printf "\nMap to chloroplast reference using bowtie2: $index"
    bowtie2 -1 $DIR/${index}_trimmed_P1.fastq -2 $DIR/${index}_trimmed_P2.fastq -x $indexedGenome -S $DIR/${index}_mapped.sam -I 200 -X 700 --very-sensitive-local --met-file $DIR/${index}_mapped.log --threads 4

    printf "\nConvert sam to bam: $index" # Converts from a human-readable format to a binary format
    samtools view -bS $DIR/${index}_mapped.sam > $DIR/${index}_mapped.bam

    printf "\nSorting reads: $index" # Sorts reads by their aligned position to reference
    samtools sort $DIR/${index}_mapped.bam > $DIR/${index}_sorted.bam

    printf "\nDetermining read depth: $index" # Outputs sequencing depth in a tab delimited format
    samtools depth $DIR/${index}_sorted.bam > $DIR/${index}_read_depth.tsv

    printf "\nFiltering by phred quality 20: $index" # Creates a pileup file; summarizes the base calls of aligned reads
    samtools mpileup -Q 20 -Agf $referenceGenome $DIR/${index}_sorted.bam > $DIR/${index}.mpilup

    printf "\nGenerating consensus genotypes: $index"
    # 
    bcftools call --ploidy 1 -c $DIR/${index}.mpilup > $DIR/${index}.vcf  #Assuming ploidy is 1 for chloroplast

    printf "\nFiltering for read depth >= 10: $index" # not sure this code works
    vcftools --minDP 10 --vcf $DIR/${index}.vcf --out $DIR/filtered --recode --recode-INFO-all

    printf "Convert vcf file back into a fasta"
    vcfutils.pl vcf2fq $DIR/filtered.recode.vcf > $DIR/${index}_assembled.fastq # -d masking for heterozygosity and low coverage
    seqtk seq -A $DIR/${index}_assembled.fastq > $DIR/${index}_assembled.fasta


    #echo "Cleaning up for now.."

        #mkdir ${output_dir}${index}
        #mv ${index}.fasta ${index}_read_depth.tsv ${output_dir}${index}
        

done