#!/bin/sh

#######################################################
## Code used to clean raw fastq files
#######################################################
## Author: Noemie Lefrancq
## Last modification: 07/03/2022
#######################################################

## WARNING: files must be fastq (paired)

## Load required modules
module load fastqc/0.11.5
module load cutadapt/1.9.1
module load pigz

## Check quality of fastq
fastqc -t 2 $( basename ${file} _R1_001.fastq)_R1_001.fastq $( basename ${file} _R1_001.fastq)_R2_001.fastq;

## Trim reads based on quality, and removed reads that are too short
cutadapt -q 30 --minimum-length=50 --pair-filter=any -o $( basename ${file} _R1_001.fastq)_1_trim.fastq -p $( basename ${file} _1.fastq.gz)_2_trim.fastq $( basename ${file} _1.fastq.gz)_R1_001.fastq $( basename ${file} _1.fastq.gz)_R2_001.fastq ;

## Check quality of the trimmed fastqc
fastqc -t 2 $( basename ${file} _R1_001.fastq)_1_trim.fastq $( basename ${file} _R1_001.fastq)_2_trim.fastq ;

## Compress output
pigz -p 2  $( basename ${file} _R1_001.fastq)_R1_001.fastq ;
pigz -p 2  $( basename ${file} _R1_001.fastq)_R2_001.fastq ;
pigz -p 2  $( basename ${file} _R1_001.fastq)_1_trim.fastq ;
pigz -p 2  $( basename ${file} _R1_001.fastq)_2_trim.fastq ;

## Move intput files to storage folder
mkdir fastq_raw
mv *_R1_001.fastq fastq_raw/
mv *_R2_001.fastq fastq_raw/

## Move fastqc reports to output folder
mkdir reports
mv $( basename ${file} _R1_001.fastq)*_fastqc.* reports/


