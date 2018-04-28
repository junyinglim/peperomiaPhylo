
#DIR=""
#cd $DIR
#ls *fastq | awk '{ print "fastq_quality_filter -p 90 -q 30" }'
#ls *.assembled.fastq | awk '{print "fastq_quality_filter -Q33 -p 90 -q 30 -i "$0" | fastq_to_fasta -Q33 -o "$0""}' | sed 's/.assembled.fastq/_loci.fasta/2' > fastqFilter.sh



#bowtie2 -x <insertfastahere?> -1 $BT2_HOME/example/reads/reads_1.fq -2 $BT2_HOME/example/reads/reads_2.fq -S eg2.sam


# Build index file
# bowtie2-build referenceGenome_KT223569.1.fasta piperchloroplast

#!/bin/bash

#rep=${num}
# loop over all samples
DIR="/Users/junyinglim/Desktop/testtest"
trimmomatic="/Users/junyinglim/Dropbox/Projects/Programs/Trimmomatic-0.36"

for s in $DIR/*_L008_R1*.fastq.gz
do
    s=${s##*/}
    index=${s%_S*}

    echo "Uncompressing index: $index"
    gunzip $DIR/${index}_*R1*.fastq.gz -k -d
    gunzip $DIR/${index}_*R2*.fastq.gz -k -d
    
    echo "Filtering adapter sequences"
    java -jar $trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 $DIR/${index}_*R1*.fastq $DIR/${index}_*R2*.fastq $DIR/${index}_trimmed_P1.fastq $DIR/${index}_trimmed_U1.fastq $DIR/${index}_trimmed_P2.fastq $DIR/${index}_trimmed_U2.fastq ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE.fa:1:30:10 SLIDINGWINDOW:10:20 MINLEN:40
    # Using a combination of the example, and Will Freyman's parameters

    #echo "Clean up intermediate files"
    #rm $DIR/${index}_*R[0-9]*.fastq

done


#     # get index from file name
#     s=${s##*/}
#     index=${s%_L001_*}

#     echo "Uncompressing index: $index"
#     gunzip data/org_sequence_data${rep}/${index}_*R1*.fastq.gz
#     gunzip data/org_sequence_data${rep}/${index}_*R2*.fastq.gz
        
#     echo "Filtering adapter sequences..."
#     java -classpath utilities/Trimmomatic-0.36/trimmomatic-0.36.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 4 -phred33

# data/org_sequence_data${rep}/${index}_*R1*.fastq

# data/org_sequence_data${rep}/${index}_*R2*.fastq

# ${index}_trimmed_P1.fq
# ${index}_trimmed_U1.fq
# ${index}_trimmed_P2.fq
# ${index}_trimmed_U2.fq ILLUMINACLIP:utilities/Fast-Plast/bin/NEB-PE.fa:1:30:10 SLIDINGWINDOW:10:20 MINLEN:40

#     # make reference guided assemblies of each target locus:
#     for f in data/reference_sequences/*.fasta  
#     do
#         # remove everything in path before slash
#         f=${f##*/} 
#         # remove file extension
#         ref=${f%.*}
#         reference="data/reference_sequences/${ref}.fasta"
#         if [ ! -d data/assemblies/$ref ]
#         then
#             mkdir data/assemblies/$ref
#         fi
#         output_dir="data/assemblies/${ref}/"
        
#         echo "Assembling index: $index"
        
#         echo "Mapping reads with bwa..."
#         bwa mem $reference -R "@RG\tID:1\tSM:1" ${index}_trimmed_P1.fq ${index}_trimmed_P2.fq > ${index}.sam
        
#         echo "Coverting SAM to BAM..."
#         samtools view -bS ${index}.sam > ${index}.bam
        
#         echo "Sorting reads..."
#         samtools sort ${index}.bam ${index}_sorted
        
#         echo "Determining read depth..."
#         samtools depth ${index}_sorted.bam > ${index}_read_depth.tsv
        
#         echo "Making BAM pileup, filtering by phred quality 20..."
#         samtools mpileup -Q 20 -Agf $reference ${index}_sorted.bam > ${index}.mpilup

#         echo "Generating consensus genotypes..."
#         bcftools view -s data/samples.txt -cg ${index}.mpilup > ${index}_temp.vcf
        
#         echo "Filtering for read depth >= 10..."
#         python src/filter_vcf_read_depth.py ${index}_temp.vcf 10

#         echo "Generating final FASTA sequence..."
#         vcfutils.pl vcf2fq ${index}_temp.vcf.filtered > ${index}.fastq
#         seqtk seq -A ${index}.fastq > ${index}.fasta
        
#         echo "Cleaning up..."
#         rm ${index}_temp.vcf ${index}_temp.vcf.filtered ${index}.mpilup
#         rm ${index}_sorted.bam ${index}.bam ${index}.sam ${index}.fastq
#         #mkdir ${output_dir}${index}
#         #mv ${index}.fasta ${index}_read_depth.tsv ${output_dir}${index}
#         mv ${index}.fasta ${index}_read_depth.tsv ${output_dir}
    
#     done

#     echo "Re-compressing index: $index"
#     rm ${index}_trimmed_*
#     gzip data/org_sequence_data${rep}/${index}_*R1*.fastq
#     gzip data/org_sequence_data${rep}/${index}_*R2*.fastq

# done

# # Illumina adapter sequences were removed and the raw sequence reads were quality filtered using Trimmomatic v0.36 (CITE). Reads were trimmed when the average phred quality score in a 10-bp sliding window was less than 20. Reads that were less than 40 bp or that did not survive the filtering process in both forward and reverse directions were excluded.
