#!/bin/sh
##file must be fastq (paired)

module load fastqc/0.11.5
module load cutadapt/1.9.1
module load pigz

## Check quality of fastq
fastqc -t 2 $( basename ${file} _R1_001.fastq)_R1_001.fastq $( basename ${file} _R1_001.fastq)_R2_001.fastq;
## TRIM READS
cutadapt -q 30 --minimum-length=50 --pair-filter=any -o $( basename ${file} _R1_001.fastq)_1_trim.fastq -p $( basename ${file} _1.fastq.gz)_2_trim.fastq $( basename ${file} _1.fastq.gz)_R1_001.fastq $( basename ${file} _1.fastq.gz)_R2_001.fastq ;
## CHECK QUALITY FASTQC TRIMMED	
fastqc -t 2 $( basename ${file} _R1_001.fastq)_1_trim.fastq $( basename ${file} _R1_001.fastq)_2_trim.fastq ;
pigz -p 2  $( basename ${file} _R1_001.fastq)_R1_001.fastq ;
pigz -p 2  $( basename ${file} _R1_001.fastq)_R2_001.fastq ;
pigz -p 2  $( basename ${file} _R1_001.fastq)_1_trim.fastq ;
pigz -p 2  $( basename ${file} _R1_001.fastq)_2_trim.fastq ;

mkdir fastq_raw
mv *_R1_001.fastq fastq_raw/
mv *_R2_001.fastq fastq_raw/

mkdir reports
mv $( basename ${file} _R1_001.fastq)*_fastqc.* reports/


