#!/bin/bash

# SCRIPT: ChIPseq-HU-SR-pipeline.sh
# AUTHOR: Mark Boltengagen
# VERSION: v16-05-2018
# OS: Linux Ubuntu 16.04 LTS
# DESCRIPTION: Transform a set of fastq.gz files to normilized bedgraph, ready for working with Bioconductor
# USED TOOLS: fastqc > trimmomatic > Bowtie2 > samtool > samtool-sort > MACS2 peaks-calling for peak calling optimized for histones

### PARAMETERS
exppath=${HOME}/Projects/007_ChIPseq_HU-SR/
Sources=${HOME}/Projects/Sources/

inpath=${exppath}raw/
outpath=${exppath}output/
mkdir -p $outpath

# How many processors are available for the computation?
threads=8

{

#####======================================================

##### RENAME RAW FILES TO CUT UNNECESSARY POSTFIXES [_S*_R*_*]

cd ${inpath}
for item in *.fastq.gz
	do
		if [[ $item = *'_S'*'_R'*'.fastq.gz' ]]; then
			mv $item ${item%*_S*_R*_*.fastq.gz}.fastq.gz
			echo "file $item was renamed into ${item%*_S*_R*_*.fastq.gz}.fastq.gz"
		fi
	done

######======================================================

##### 1. RAW READS COUNTER
##### Count the number of total and unique reads, the most abundat sequences which might be a possible contamination
##### Reference: Bardet et al., Nature Protocols, 2011

cd ${inpath}
printf "$(date)\tSTARTED STATISTICS OF RAW READS\n\n"
printf "Name\tTotalReads\tUnique\tUnique%%\tMostAbundantSeq\tScore\tScore%%\n"
for item in *.fastq.gz
	do
		gunzip -c $item | awk -v var="$item" '((NR-2)%4==0) {read=$1;total++;count[read]++}END{for(read in count) {if(!max||count[read]>max){max=count[read];maxRead=read};if(count[read]==1){unique++}};print var, total, unique, unique*100/total,maxRead,count[maxRead], count[maxRead]*100/total}'
	done
printf "=%.0s" {1..70}; printf "\n\n"

#####======================================================

##### 2. FASTQC FASTQ QUALITY TEST

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

##### 3. TRIMMOMATIC

trimmomatic_path="/home/suvar/Softs/Trimmomatic-0.36/trimmomatic-0.36.jar"

cd ${inpath}
mkdir ${outpath}trim
printf "$(date)\tSTARTED TRIMMING OF ADAPTERS\n\n"
for item in *.fastq.gz
	do
		java -jar $trimmomatic_path SE -phred33 $item ${outpath}trim/${item%.fastq.gz}-trim.fastq.gz HEADCROP:12 SLIDINGWINDOW:6:20 CROP:40 MINLEN:40
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
		fastqc -t $threads -o ${outpath}fastQCtrim/ $item
		printf "\n"
	done
printf "$(date)\tFINISHED FASTQ OF TRIMMED FASTQ\n\n"
printf "=%.0s" {1..80}; printf "\n\n"

#####======================================================

##### 5. TRIMMED READS COUNTER
## Count the number of total and unique reads, the most abundat sequences which might be a possible contamination
## Reference: Bardet et al., Nature Protocols, 2011

cd ${outpath}trim
printf "$(date)\tSTARTED STATISTICS OF TRIMMED READS\n\n"
printf "Name\tTotalReads\tUnique\tUnique%%\tMostAbundantSeq\tScore\tScore%%\n"
for item in *.fastq.gz
	do
		gunzip -c $item | awk -v var="$item" '((NR-2)%4==0) {read=$1;total++;count[read]++}END{for(read in count) {if(!max||count[read]>max){max=count[read];maxRead=read};if(count[read]==1){unique++}};print var, total, unique, unique*100/total,maxRead,count[maxRead], count[maxRead]*100/total}'
	done
printf "\n$(date)\tCOMPLETED STATISTICS OF TRIMMED READS\n\n"
printf "=%.0s" {1..80}; printf "\n\n"

#####======================================================

##### 6. BOWTIE2 for SacCer3
indexSacCer3="${Sources}S288C_SacCer3_index/S288C_SacCer3_index"

cd ${outpath}trim
mkdir ${outpath}bam
printf "$(date)\tBOWTIE2 FOR SACCER3 STARTS\n"

for item in *.fastq.gz
	do
		printf "\n$(date)\tStarted alignment of ${item}\n"
		bowtie2 -p $threads -x ${indexSacCer3} -U $item | samtools view -@ $threads -bS - > ${outpath}bam/${item%-trim.fastq.gz}.bam
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

#####======================================================

#### 7. BOWTIE2 for Spikes (Pombe)

indexPombe="${Sources}pombe_index/pombe_index"

cd ${outpath}trim
mkdir ${outpath}bam-spike
printf "$(date)\tBOWTIE2 FOR SPIKES STARTS\n"

for item in *.fastq.gz
	do
		printf "\n$(date)\tStarted alignment of ${item}\n"
		bowtie2 -p $threads -x ${indexPombe} -U $item | samtools view -@ $threads -bS - > ${outpath}bam-spike/${item%-trim.fastq.gz}.bam
		printf "$(date)\tCompleted alignment of ${item}\n"
	done
printf "=%.0s" {1..80}; printf "\n\n"

cd ${outpath}bam-spike
mkdir ${outpath}fbam-spike
printf "$(date)\tSAMTOOLS FILTERING FOR SPIKES STARTS\n"

for item in *.bam
	do
		printf "\n$(date)\tStarted bam filtering of ${item}\n"
		samtools view -@ $threads -F 4 -u -q 20 $item | samtools sort -@ $threads - | samtools rmdup -s - ${outpath}fbam-spike/${item%.bam}.bam
		printf "$(date)\tCompleted bam filtering of ${item}\n"
	done
printf "=%.0s" {1..80}; printf "\n\n"

#####======================================================

##### MAPPED READS STATISTICS
##### Adopted from Bardet et al., Nature Protocols, 2011

cd ${outpath}bam

printf "\n\nName \tRAW \tTotalReads \tTotalReads%% \tUnique \tUnique%% \tMaxCoordinates \tMax%%\n"

for item in *.bam
	do
		# Number of raw reads
		raw=$(samtools view $item | wc -l)
		# Number of raw, unique and most repeated reads (bedtools)
		bedtools bamtobed -i $item | awk -v RAW=$raw -v name=$item '{coordinates=$1":"$2"-"$3;total++;count[coordinates]++}END{for(coordinates in count){if(!max||count[coordinates]>max){max=count[coordinates];maxCoor=coordinates};if(count[coordinates]==1){unique++}};print name,RAW,total,total*100/RAW,unique,unique*100/total,maxCoor,count[maxCoor]*100/total}'

		# Total and top 10 of non-mapped reads
		samtools view -f 0x0004 $item | awk '{reads=$10;total++;count[reads]++}END{print "Total_non-mapped_reads",total;for(read in count){print read,count[read]+0}}' | sort -k2,2nr |head -11
	done

##### MAPPED READS STATISTICS FOR SPIKES
##### Adopted from Bardet et al., Nature Protocols, 2011

cd ${outpath}bam-spike

printf "\n\nName \tRAW \tTotalReads \tTotalReads%% \tUnique \tUnique%% \tMaxCoordinates \tMax%%"

for item in *.bam
	do
		# Number of raw reads
		raw=$(samtools view $item | wc -l)
		# Number of raw, unique and most repeated reads
		bedtools bamtobed -i $item | awk -v RAW=$raw -v name=$item '{coordinates=$1":"$2"-"$3;total++;count[coordinates]++}END{for(coordinates in count){if(!max||count[coordinates]>max){max=count[coordinates];maxCoor=coordinates};if(count[coordinates]==1){unique++}};print name,RAW,total,total*100/RAW,unique,unique*100/total,maxCoor,count[maxCoor]*100/total}'

		# Total and top 10 of non-mapped reads
		samtools view -f 0x0004 $item | awk '{reads=$10;total++;count[reads]++}END{print "Total_non-mapped_reads",total;for(read in count){print read,count[read]+0}}' | sort -k2,2nr |head -11
	done

#####======================================================

##### FINDING PAIRS IP_N ~ INP_N

cd ${outpath}fbam
mkdir ${outpath}macs2

for ip in *IP.bam
	do	
		pre=$(echo $ip | awk -F "_" '{print $2}')
		exp=$(echo $ip | awk -F "_" '{print $3}')
		time=$(echo $ip | awk -F "_" '{print $4}')
		input=$(find -name *${pre}_${exp}_${time}_IN.bam)
		echo ip sample: $ip
		echo input sample: $input

##### MACS2 PEAK CALLING

		printf "$(date)\tMacs2 callpeak for input file - $input and IP file - $ip starts\n"
		macs2 callpeak --nomodel --keep-dup all --shift 37 --extsize 73 --broad -f BAM -g12157105 -q 0.05 -t $ip -c $input --outdir ${outpath}macs2 -n ${ip%.bam} -B --SPMR
		printf "$(date)\tMacs2 callpeak for input file - $input and IP file - $ip is done\n\n"
	done
printf "=%.0s" {1..80}; printf "\n\n"

#####======================================================

##### ADDING TRACK LINE INTO BEDGRAPH FILES

cd ${outpath}macs2

printf "$(date)\tBedGraph track preparation starts\n"

for item in *treat_pileup.bdg
	do
		printf "$(date)\tStarted bedGraph track preparation of ${item}\n"
		sed -i '1 i\track type=bedGraph name='"'"$item"'"' visibility=FULL color=0,0,255 graphType=bar' $item
		printf "$(date)\tCompleted bedGraph track preparation of ${item}\n\n"
	done
printf "=%.0s" {1..80}; printf "\n\n"

#####======================================================

##### BDG SORTING

cd ${outpath}macs2
mkdir ${outpath}bdg
## Remove header if exist and sort bdg file
for item in *treat_pileup.bdg
	do
		sed -e '1{/^track type=bedGraph/d}' $item | LC_COLLATE=C sort -k1,1 -k2,2n | gzip > ${outpath}bdg/${item}.gz 
	done

#####=======================================================

##### COMPRESS MACS2 bedGraphs
printf "$(date)\MACS2 result compression starts\n"
cd ${outpath}macs2
for item in *
	do
		printf "$(date)\tbedGraph compression of ${item}\n"
		gzip ${item}
		printf "$(date)\tCompleted bedGraph compression of ${item}\n\n"
	done
printf "=%.0s" {1..80}; printf "\n\n"

#####=======================================================

##### STATISTICS COLLECTION

trim_dir=${outpath}trim/
mapped_dir=${outpath}bam/
duprmv_dir=${outpath}fbam/
mapped_dir_spike=${outpath}bam-spike/
duprmv_dir_spike=${outpath}fbam-spike/

printf "Name\tTotal_reads\tTrimmed_reads\t%%Trimmed\tMapped_reads\t%%_Mapped\tFiltered_reads\t%%_Filtered\tMapped_reads_spikes\t%%_Mapped_spikes\tFiltered_reads_spikes\t%%_Filtered_spikes\n" >> ${outpath}statistics.tsv
cd ${inpath}
for item in *.fastq.gz
	do
		prefix=$(echo $item | awk -F "_" '{print $1}')

		total_reads=$(gunzip -c $item | awk '(NR-2)%4 ==0 {total++}END{print total}')
		trim_reads=$(gunzip -c ${trim_dir}${prefix}_*.fastq.gz | awk '(NR-2)%4 ==0 {total++}END{print total}')
		rate_trim=$(echo "scale=2; $trim_reads/$total_reads" | bc)
		mapped_reads=$(samtools view -c -F 4 ${mapped_dir}${prefix}_*.bam)
		rate_mapped=$(echo "scale=2; $mapped_reads/$trim_reads" | bc)
		duprmv_reads=$(samtools view -c ${duprmv_dir}${prefix}_*.bam)
		rate_filtered=$(echo "scale=2; $duprmv_reads/$mapped_reads" | bc)
		mapped_reads_spike=$(samtools view -c -F 4 ${mapped_dir_spike}${prefix}_*.bam)
		rate_mapped_spike=$(echo "scale=2; $mapped_reads_spike/$trim_reads" | bc)
		duprmv_reads_spike=$(samtools view -c ${duprmv_dir_spike}${prefix}_*.bam)
		rate_filtered_spike=$(echo "scale=2; $duprmv_reads_spike/$mapped_reads_spike" | bc)

		printf ${item}'\t'${total_reads}'\t'${trim_reads}'\t'${rate_trim}'\t'${mapped_reads}'\t'${rate_mapped}"\t"${duprmv_reads}'\t'${rate_filtered}"\t"${mapped_reads_spike}'\t'${rate_mapped_spike}"\t"${duprmv_reads_spike}'\t'${rate_filtered_spike}'\n' >> ${outpath}statistics.tsv
	done
printf "=%.0s" {1..70}; printf "\n\n"

##### BDG to BigWig
cd ${outpath}bdg
# unzip bdg files, keep .gz archives
for item in *.bdg.gz; do gunzip -k $item; done

mkdir ${outpath}bw

##### BedGraph to BigWig conversion
for item in *.bdg
	do
		printf "$(date)\tBedGraph to BigWig conversion of ${item}\n"
        /home/suvar/Softs/UCSC/bedGraphToBigWig ${item} ${Sources}sacCer3.chrom.sizes ${outpath}bw/${item%_treat_pileup.bdg}.bw
        printf "$(date)\tCompleted BedGraph to BigWig conversion of ${item}\n\n"

	done
printf "=%.0s" {1..70}; printf "\n\n"

# remove unzipped bdg files
for item in *.bdg; do rm $item; done

printf "$(date)\tPIPELINE IS COMPLETED\n\n"
} 2>&1 | tee ${outpath}ChIPseq-HU-SR-pipeline.log
