#!/bin/bash

# SCRIPT: SNP_workflow.sh
# AUTHOR: Mark Boltengagen
# VERSION: v29-11-2018
# OS: Linux Ubuntu 18.04 LTS
# DESCRIPTION: SNP analysis
# USED TOOLS: bcftools

### PARAMETERS
exppath=${HOME}/Projects/011_SNP/output/
Sources=${HOME}/Projects/Sources/

# How many processors are available for the computation?
threads=4

{

##### 1. Variation call
index="${Sources}SacCer3_index_samtools/S288C_SacCer3_genome.fasta"

cd ${exppath}fbam
mkdir ${exppath}flt_vcf
printf "$(date)\tVARIATION CALL STARTED\n"

for item in *.bam
	do
		printf "\n$(date)\tStarted filtering of ${item}\n"
		samtools mpileup -uf $index $item | \
        bcftools call -mv -P 1 | \
        bcftools filter -Oz -s LowQual -e '%QUAL<20 || DP>100'  - -o ${exppath}flt_vcf/${item%.bam}.var.flt.vcf.gz
        printf "$(date)\tCompleted filtering of ${item}\n"
	done

##### 2. Compaire two vcf files
cd ${exppath}flt_vcf
mkdir ${exppath}dif_vcf
printf "$(date)\tCOMPAIR VCF FILES STARTED\n"

## make index of created vcf files for downstream analysis
for item in *.var.flt.vcf.gz
	do
		bcftools index $item
	done

printf "$(date)\tStarted compair of files\n"
bcftools isec -c indels $(ls *.var.flt.vcf.gz) -p ${exppath}dif_vcf
printf "$(date)\tCompleted compair of files\n"

printf "$(date)\tPIPELINE IS COMPLETE\n\n"
} 2>&1 | tee ${outpath}SNP_pipeline.log
