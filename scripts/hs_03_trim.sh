#!/bin/bash
# this script loads trimmomatic, trims the reads, and moves the output to
# data/fastq_trim
log_file="$(date "+%m%d%Y-%H%M%S").log"

echo "Loading Trimmomatic"
module load java

cd ../data/fastq
echo "Working in" $PWD
# iterating through the first file of each paired end set of .fastq.gz files
for infile in *_1.fastq.gz
do
base=$(basename ${infile} _1.fastq.gz)
echo "Trimming "$base
# change number of threads as necessary for computing hardware
java -jar ~/modules/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 19 ${infile} ${base}_2.fastq.gz ${base}_1.trim.fastq.gz ${base}_1.orphan.fastq.gz ${base}_2.trim.fastq.gz ${base}_2.orphan.fastq.gz SLIDINGWINDOW:4:20 MINLEN:20 ILLUMINACLIP:~/modules/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:40:15 > $log_file 
done

echo "Moving trimmed fastq files"
mkdir ../fastq_trim
mv *.trim.fastq.gz ../fastq_trim/
