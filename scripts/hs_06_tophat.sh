#!/bin/bash
# this script runs tophat on the trimmed fastq files

log_file="$(date "+%m%d%Y-%H%M%S").log"

# create variables for reference genome and annotation file locations relative to directory
# where the mapping will be done
reference=../../ref/dm6
annotation=../../ref/Drosophila_melanogaster.BDGP6.95.gtf

# make a folder for mapping results
mkdir ../data/tophat_map
cd ../data/tophat_map

# do the mapping
for f in *_1.trim.fastq.gz
do
sample=$(basename ${f} _1.trim.fastq.gz)
echo "Running tophat on" $sample
# change number of threads as appropriate for computing hardware
tophat2 -p 16 -G ${annotation} -o ${sample} ${reference} ../fastq_trim/${sample}_1.trim.fastq.gz ../fastq_trim/${sample}_2.trim.fastq.gz
done > $log_file
