#!/bin/bash

# SCRIPT: NGS_workflow.sh
# AUTHOR: Mark Boltengagen
# VERSION: v28-11-2018
# OS: Linux Ubuntu 18.04 LTS
# DESCRIPTION: SNP analysis
# USED TOOLS: fastqc > Bowtie2 > samtool > bcftools

### PARAMETERS
exppath=${HOME}/Projects/011_SNP/
Sources=${HOME}/Projects/Sources/

inpath=${exppath}raw/
outpath=${exppath}output/
mkdir -p $outpath

# How many processors are available for the computation?
threads=16

{

######======================================================

##### 1. FASTQC FASTQ QUALITY TEST

cd ${inpath}
mkdir ${outpath}/fastQC
printf "$(date)\tSTARTED FASTQ OF RAW FASTQ\n\n"
for item in *.fastq.gz
	do
		fastqc -t $threads -o ${outpath}fastQC/ $item
		printf "\n"
	done
printf "$(date)\tFINISHED FASTQ OF RAW FASTQ\n\n"
printf "=%.0s" {1..80}; printf "\n\n"

#####======================================================

##### 2. BOWTIE2 for SacCer3 -PE
indexSacCer3="${Sources}S288C_SacCer3_index/S288C_SacCer3_index"

cd ${inpath}
mkdir ${outpath}bam
printf "$(date)\tBOWTIE2 FOR SACCER3 STARTS\n"

for item in *_1.fastq.gz
	do
		printf "\n$(date)\tStarted alignment of ${item}\n"
		bowtie2 -p $threads -x ${indexSacCer3} -1 $item -2 ${item%1.fastq.gz}2.fastq.gz | samtools view -@ $threads -bS - > ${outpath}bam/${item%_1.fastq.gz}.bam
		printf "$(date)\tCompleted alignment of ${item}\n"
	done
printf "=%.0s" {1..80}; printf "\n\n"

cd ${outpath}bam
mkdir ${outpath}fbam
printf "$(date)\tSAMTOOLS FILTERING FOR SACCER3 STARTS\n"

for item in *.bam
	do
		printf "\n$(date)\tStarted bam filtering of ${item}\n"
		samtools view -@ $threads -F 4 -u -q 20 $item | samtools sort -@ $threads - | samtools rmdup -s - ${outpath}fbam/${item%.bam}.bam
		printf "$(date)\tCompleted bam filtering of ${item}\n"
	done
printf "=%.0s" {1..80}; printf "\n\n"

printf "$(date)\tPIPELINE IS COMPLETED\n\n"

} 2>&1 | tee ${outpath}NGS_workflow.log
