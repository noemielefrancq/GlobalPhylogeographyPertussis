#!/bin/sh

#######################################################
## Script used to map reads and perform
##            variant calling 
#######################################################
### Author: Noemie Lefrancq
### Last modification: 07/03/2022
#######################################################

## Load required modules
module load bwa/0.7.7
module load java/1.8.0
module load picard-tools/2.8.1
module load samtools/1.9


#################################################### 
## First step: Mapping with BWA
#################################################### 
# Create pe.sam files out of trimmed fastq

# Unzip input files
pigz -p 4 -d $( basename ${file} _1_trim.fastq.gz)_1_trim.fastq.gz ;
pigz -p 4 -d $( basename ${file} _1_trim.fastq.gz)_2_trim.fastq.gz ;

# Run bwa mem
R1=$( basename ${file} _1_trim.fastq.gz)""_1_trim.fastq
R2=$( basename ${file} _1_trim.fastq.gz)""_2_trim.fastq
bwa mem -t 4 TohamaNC_0029292.fasta $R1 $R2 > $( basename ${file} _1_trim.fastq.gz)_trim.pe.sam

# Zip input files again
pigz -p 4 $R1 ;
pigz -p 4 $R2 ;

#################################################### 
## Second step: Variant calling with GATK
##################################################### 
mkdir tmp_$( basename ${file} _1_trim.fastq.gz)
export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp_$( basename ${file} _1_trim.fastq.gz)

# Sort sam file
/pasteur/homes/nlefranc/software/gatk-4.1.4.0/./gatk SortSam -I $( basename ${file} _1_trim.fastq.gz)_trim.pe.sam -O $( basename ${file} _1_trim.fastq.gz)""_trim_sorted.bam -SO coordinate
rm $( basename ${file} _1_trim.fastq.gz)_trim.pe.sam

# Mark duplicates, add read group and create a bam index
gatk-4.1.4.0/./gatk MarkDuplicates -I $( basename ${file} _1_trim.fastq.gz)""_trim_sorted.bam -O $( basename ${file} _1_trim.fastq.gz)""_trim_dedup.bam -M $( basename ${file} _1_trim.fastq.gz)""metrics.txt 
picard AddOrReplaceReadGroups I=$( basename ${file} _1_trim.fastq.gz)""_trim_dedup.bam O=$( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
gatk-4.1.4.0/./gatk BuildBamIndex -I $( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bam

# Run haplotypecaller
gatk-4.1.4.0/./gatk HaplotypeCaller -R TohamaNC_0029292.fasta -I $( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bam -ploidy 1 -ERC GVCF --annotation StrandBiasBySample -O $( basename ${file} _1_trim.fastq.gz)_GATK.vcf

# Remove tmp directory
rm -r tmp_$( basename ${file} _1_trim.fastq.gz)

##################################################### 
## Clean: move output files in folers
##################################################### 
mkdir vcf
mkdir bam

mv $( basename ${file} _1_trim.fastq.gz)""_GATK.vcf vcf/
mv $( basename ${file} _1_trim.fastq.gz)""_trim_sorted.bam bam/

rm $( basename ${file} _1_trim.fastq.gz)""_trim_dedup.bam 
rm $( basename ${file} _1_trim.fastq.gz)""metrics.txt 
rm $( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bam 
rm $( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bai
rm $( basename ${file} _1_trim.fastq.gz)""_GATK.vcf.idx;
