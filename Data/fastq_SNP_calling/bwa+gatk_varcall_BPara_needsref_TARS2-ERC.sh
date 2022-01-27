#!/bin/sh

#########################################################
## BWA mapping + Variant calling with GATK
## from Nadia's script
###LATEST###17Oct2019
module load bwa/0.7.7
module load java/1.8.0
module load picard-tools/2.8.1
module load samtools/1.9

#################################################### Mapping with BWA

#Create pe.sam files out of trimmed fastq
pigz -p 4 -d $( basename ${file} _1_trim.fastq.gz)_1_trim.fastq.gz ;
pigz -p 4 -d $( basename ${file} _1_trim.fastq.gz)_2_trim.fastq.gz ;
R1=$( basename ${file} _1_trim.fastq.gz)""_1_trim.fastq
R2=$( basename ${file} _1_trim.fastq.gz)""_2_trim.fastq
bwa mem -t 4 /pasteur/scratch/users/nlefranc/refseq/TohamaNC_0029292.fasta $R1 $R2 > $( basename ${file} _1_trim.fastq.gz)_trim.pe.sam
pigz -p 4 $R1 ;
pigz -p 4 $R2 ;

##################################################### Variant calling
mkdir tmp_$( basename ${file} _1_trim.fastq.gz)
export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp_$( basename ${file} _1_trim.fastq.gz)

#Sort sam file
/pasteur/homes/nlefranc/software/gatk-4.1.4.0/./gatk SortSam -I $( basename ${file} _1_trim.fastq.gz)_trim.pe.sam -O $( basename ${file} _1_trim.fastq.gz)""_trim_sorted.bam -SO coordinate
rm $( basename ${file} _1_trim.fastq.gz)_trim.pe.sam

#Mark duplicates, add read group and create a bam index
/pasteur/homes/nlefranc/software/gatk-4.1.4.0/./gatk MarkDuplicates -I $( basename ${file} _1_trim.fastq.gz)""_trim_sorted.bam -O $( basename ${file} _1_trim.fastq.gz)""_trim_dedup.bam -M $( basename ${file} _1_trim.fastq.gz)""metrics.txt 
picard AddOrReplaceReadGroups I=$( basename ${file} _1_trim.fastq.gz)""_trim_dedup.bam O=$( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
/pasteur/homes/nlefranc/software/gatk-4.1.4.0/./gatk BuildBamIndex -I $( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bam

/pasteur/homes/nlefranc/software/gatk-4.1.4.0/./gatk HaplotypeCaller -R /pasteur/scratch/users/nlefranc/refseq/TohamaNC_0029292.fasta -I $( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bam -ploidy 1 -ERC GVCF --annotation StrandBiasBySample -O $( basename ${file} _1_trim.fastq.gz)_GATK.vcf

rm -r tmp_$( basename ${file} _1_trim.fastq.gz)
##################################################### Compute coverage
## compute distribution of coverage
samtools stats --coverage 5,1000,1  $( basename ${file} _1_trim.fastq.gz)""_trim_sorted.bam | grep ^COV | cut -f 2- >  $( basename ${file} _1_trim.fastq.gz)""_trim_coverage_distribution.txt

## compute coverage at each position
samtools depth -aa $( basename ${file} _1_trim.fastq.gz)""_trim_sorted.bam > $( basename ${file} _1_trim.fastq.gz)""_trim_coverage_each_position_tmp.txt ## computes a long file containing various statistics
cut -f 3 $( basename ${file} _1_trim.fastq.gz)""_trim_coverage_each_position_tmp.txt > $( basename ${file} _1_trim.fastq.gz)""_trim_coverage_each_position.txt ## keeps only the coverage at each position

mkdir coverage
rm $( basename ${file} _1_trim.fastq.gz)""_trim_coverage_each_position_tmp.txt
mv $( basename ${file} _1_trim.fastq.gz)""_trim_coverage_distribution.txt $( basename ${file} _1_trim.fastq.gz)""_trim_coverage_each_position.txt coverage

##################################################### Clean: put files in folers
mkdir vcf
mkdir bam

mv $( basename ${file} _1_trim.fastq.gz)""_GATK.vcf vcf/
mv $( basename ${file} _1_trim.fastq.gz)""_trim_sorted.bam bam/

rm $( basename ${file} _1_trim.fastq.gz)""_trim_dedup.bam 
rm $( basename ${file} _1_trim.fastq.gz)""metrics.txt 
rm $( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bam 
rm $( basename ${file} _1_trim.fastq.gz)""_trim_dedup2.bai
rm $( basename ${file} _1_trim.fastq.gz)""_GATK.vcf.idx;
