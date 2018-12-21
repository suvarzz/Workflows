#!/bin/bash

# SCRIPT: RNAseq-pipeline.sh
# AUTHOR: Mark Boltengagen
# VERSION: v14-05-2018
# OS: Linux Ubuntu 16.04 LTS
# DESCRIPTION: RNA-seq pipeline

##### PARAMETERS
inpath=${HOME}/Projects/003_RNAseq/raw/
outpath=${HOME}/Projects/003_RNAseq/output/
sourcepath=${HOME}/Projects/Sources/
mkdir -p $outpath

# How many processors are available for the computation?
threads=4

{

##### RENAME RAW FILES TO CUT UNNECESSARY POSTFIXES [_S*_R*_*]

cd ${inpath}
for item in *.fastq.gz
	do
		if [[ $item = *'_S'*'_R'*'.fastq.gz' ]]; then
			mv $item ${item%*_S*_R*_*.fastq.gz}.fastq.gz
			echo "file $item was renamed into ${item%*_S*_R*_*.fastq.gz}.fastq.gz"
		fi
	done

#####======================================================

##### 2. FASTQC FASTQ QUALITY TEST

cd ${inpath}
mkdir ${outpath}/fastQC
printf "$(date)\tSTARTED FASTQ OF RAW FASTQ\n\n"
for item in *.fastq.gz
	do
		fastqc -o ${outpath}fastQC/ $item
		printf "\n"
	done
printf "$(date)\tFINISHED FASTQ OF RAW FASTQ\n\n"
printf "=%.0s" {1..80}; printf "\n\n"

#####======================================================

##### 3. TRIMMOMATIC

trimmomatic_path="/home/suvar/Softs/Trimmomatic-0.36/trimmomatic-0.36.jar"

cd ${inpath}RAW
mkdir ${outpath}trim
printf "$(date)\tSTARTED TRIMMING OF ADAPTERS\n\n"
for item in *.fastq.gz
	do
		java -jar $trimmomatic_path SE -phred33 $item ${outpath}trim/${item%.fastq.gz}-trim.fastq.gz  HEADCROP:12 SLIDINGWINDOW:6:20 CROP:50 MINLEN:50
		printf "\n"
	done
printf "$(date)\tFINISHED TRIMMING OF ADAPTERS\n\n"
printf "=%.0s" {1..80}; printf "\n\n"

#####======================================================

##### 4. FASTQC OF TRIMMED FILES

cd ${outpath}trim
mkdir ${outpath}fastQCtrim
printf "$(date)\tSTARTED FASTQ OF TRIMMED FASTQ\n\n"
for item in *.fastq.gz
	do
		fastqc -o ${outpath}fastQCtrim/ $item
		printf "\n"
	done
printf "$(date)\tFINISHED FASTQ OF TRIMMED FASTQ\n\n"
printf "=%.0s" {1..80}; printf "\n\n"

#####======================================================

##### Salmon alignment and count

transcripts_index=${sourcepath}/SacCer3_transcripts_index/transcripts_index/

mkdir ${outpath}salmon
cd ${outpath}trim

for item in *.fastq.gz
	do
		printf "\n$(date)\tStarted Salmon alignment and count of ${item}\n"
		salmon quant -i $transcripts_index -l A -r $item -o ${outpath}salmon/${item%.fastq.gz}
		printf "\n$(date)\tCompleted Salmon alignment and count of ${item}\n"
	done

printf "$(date)\tPIPELINE IS COMPLETE\n\n"
} 2>&1 | tee ${outpath}RNAseq_pipeline.log
