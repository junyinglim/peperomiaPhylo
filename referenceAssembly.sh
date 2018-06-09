#!/bin/bash

# Dependencies:
#   samtools, bowtie2, trimmomatic, gunzip
#   bcftools (tools for variant calling)

method="bwa" # bwa or bowtie2
genome="ribosome" # chloroplast or ribosome
hostname=$(hostname)

# Directories for the server (space between if and square bracket is necessary)
if [ "$hostname" = "voldemort" ]; then
    DIR="/home/junyinglim/peperomiaFiles"
    trimmomatic="/home/junyinglim/Trimmomatic-0.36"

    if [ "$genome" = "chloroplast" ]; then
        indexedGenome="/home/junyinglim/referenceGenome/piperchloroplast"
        referenceGenome="/home/junyinglim/referenceGenome/referenceGenome_KT223569.1.fasta"
    elif [ "$genome" = "ribosome" ]; then
        referenceGenome="/home/junyinglim/referenceGenome/JYL-31_rDNA.fasta"
    fi    
    
# Directories for local
elif [ "$hostname" = "Juns-MacBook-Pro.local" ]; then
    DIR="/Users/junyinglim/Desktop/testtest"
    trimmomatic="/Users/junyinglim/Dropbox/Projects/Programs/Trimmomatic-0.36"
    
    if [ "$genome" = "chloroplast" ]; then
        indexedGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/piperchloroplast"
        referenceGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/referenceGenome_KT223569.1.fasta"    
    elif [ "$genome" = "ribosome" ]; then
        referenceGenome="/Users/junyinglim/Dropbox/Projects/2015/Peperomia/peperomiaPhylo/referenceGenome/JYL-31_rDNA.fasta"
    fi

fi


# Index reference genomes
# if [ "$method" = "bwa" ]; then
#     bwa index $referenceGenome
# fi

for s in $DIR/*_L008_R1*.fastq.gz
do
    s=${s##*/}
    index=${s%_S*}

    # # Trimmming sequences
    # printf "\nUncompressing index: $index"
    # gunzip $DIR/${index}_*R1*.fastq.gz -k -d
    # gunzip $DIR/${index}_*R2*.fastq.gz -k -d
    
    # printf "\nFiltering adapter sequences: $index"
    # java -jar $trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 $DIR/${index}_*R1*.fastq $DIR/${index}_*R2*.fastq $DIR/${index}_trimmed_P1.fastq $DIR/${index}_trimmed_U1.fastq $DIR/${index}_trimmed_P2.fastq $DIR/${index}_trimmed_U2.fastq ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE.fa:1:30:10 SLIDINGWINDOW:10:20 MINLEN:40
    

    # USING BWA ==================
    if [ "$method" = "bwa" ]; then
        printf "\nMap to reference sequence using BWA: $index\n"
        bwa mem $referenceGenome $DIR/${index}_trimmed_P1.fastq $DIR/${index}_trimmed_P2.fastq > $DIR/${index}_bwa_${genome}_mapped.sam -t 4

        printf "\nSort reads by their aligned poisition to reference\n"
        samtools sort $DIR/${index}_bwa_${genome}_mapped.sam > $DIR/${index}_bwa_${genome}_sorted.sam

        printf "\nConvert to pile up format\n"
        samtools mpileup -Q 33 -Agf $referenceGenome $DIR/${index}_bwa_${genome}_sorted.sam > $DIR/${index}_bwa_${genome}.mpilup

        printf "\nGenerate consensus genotype: $index\n"
        bcftools call -c $DIR/${index}_bwa_${genome}.mpilup > $DIR/${index}_bwa_${genome}.vcf

        printf "\nConvert vcf file back to a fasta\n"
        vcfutils.pl vcf2fq $DIR/${index}_bwa_${genome}.vcf > $DIR/${index}_bwa_${genome}_assembled.fastq
        seqtk seq -A $DIR/${index}_bwa_${genome}_assembled.fastq > $DIR/${index}_bwa_${genome}_assembled.fasta
    fi

    # USING BOWTIE ==================
    if [ "$method" = "bowtie" ]; then
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
    fi

done

# Clean up workspace
# printf "Cleaning up intermediate files ..."
# rm *_assembled.fastq
# rm *_mapped.sam
# rm *_sorted.bam
# rm *_sorted.sam
# rm *.vcf
# rm *.mpilup