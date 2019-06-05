#!/bin/bash
# this script runs fastqc on trimmed data

log_file="$(date "+%m%d%Y-%H%M%S").log"

echo "Loading fastqc"
module load fastqc

cd ../data/fastq_trim
echo "Entering" $PWD

fastqc *.trim.fastq.gz > $log_file

echo "Moving fastqc output"
mkdir ../fastqc_trimmed
mv *.zip ../fastqc_trimmed/
mv *.html ../fastqc_trimmed/