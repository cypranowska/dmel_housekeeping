#!/bin/bash
# this script compresses each .fastq file with gzip, loads fastqc as a module, 
# then runs fastqc on each .fastq.gz file

log_file="$(date "+%m%d%Y-%H%M%S").log"

cd ../data/fastq
echo "Entering" $PWD

for f in *.fastq
do
	echo "Zipping" $f
	gzip $f
done

module load fastqc
echo "Loading fastqc"
fastqc *.fastqc.gz > $log_file

mkdir ../fastqc_raw
echo "Moving fastqc output"
mv *.zip ../fastqc_raw/
mv *.html ../fastqc_raw/
