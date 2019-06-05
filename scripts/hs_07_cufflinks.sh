#!/bin/bash
# this script runs cufflinks to assemble transcriptomes for each sample

log_file="$(date "+%m%d%Y-%H%M%S").log"

# make a folder for assembled transcripts
mkdir ../data/cufflinks_out
cd ../data/cufflinks_out

# do the mapping
module load boost bowtie2 cufflinks

for d in ../tophat_map/*/
do
sample=$(basename "${d}") 
mkdir "${sample}"
cd "${sample}"
echo "Running cufflinks on" $sample
cufflinks -v ../../tophat_map/${sample}/accepted_hits.bam
cd .. 
done > $log_file
