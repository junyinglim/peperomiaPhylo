#!/bin/bash

# Dependencies: samtools, bowtie2, trimmomatic, gunzip, bcftools


DIR="/Users/junyinglim/Desktop/testtest"
trimmomatic="/Users/junyinglim/Dropbox/Projects/Programs/Trimmomatic-0.36"
indexedGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/piperchloroplast"
referenceGenome=""

for s in $DIR/*_L008_R1*.fastq.gz
do
    s=${s##*/}
    index=${s%_S*}

    echo "Uncompressing index: $index"
    gunzip $DIR/${index}_*R1*.fastq.gz -k -d
    gunzip $DIR/${index}_*R2*.fastq.gz -k -d
    
    echo "Filtering adapter sequences: $index"
    java -jar $trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 $DIR/${index}_*R1*.fastq $DIR/${index}_*R2*.fastq $DIR/${index}_trimmed_P1.fastq $DIR/${index}_trimmed_U1.fastq $DIR/${index}_trimmed_P2.fastq $DIR/${index}_trimmed_U2.fastq ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE.fa:1:30:10 SLIDINGWINDOW:10:20 MINLEN:40

    echo "Map to chloroplast reference using bowtie2: $index"
    bowtie2 -1 $DIR/${index}_trimmed_P1.fastq -2 $DIR/${index}_trimmed_P2.fastq -x $indexedGenome -S $DIR/${index}_mapped.sam

    echo "Convert sam to bam: $index"
    # not sure why this is necessary
    samtools view -bS $DIR/${index}_mapped.sam > $DIR/${index}_mapped.bam

    echo "Sorting reads: $index"
    samtools sort $DIR/${index}_mapped.bam > $DIR/${index}_sorted.bam

    echo "Determining read depth: $index"
    samtools depth $DIR/${index}_sorted.bam > $DIR/${index}_read_depth.tsv

    # echo "Filtering by phred quality 20: $index"
    # samtools mpileup -Q 20 -Agf referenceGenome/referenceGenome_KT223569.1.fasta ${index}_sorted.bam > ${index}.mpilup

    # echo "Generating consensus genotypes..."
    # bcftools view -s data/samples.txt

    # echo "Filtering for read depth >= 10..."

    #echo "Cleaning up for now.."

        #mkdir ${output_dir}${index}
        #mv ${index}.fasta ${index}_read_depth.tsv ${output_dir}${index}
        

done