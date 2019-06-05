#!/bin/bash
# this script quantifies the abundances of each transcript with cuffquant

log_file="$(date "+%m%d%Y-%H%M%S").log"

#create a variable with paths to cuffmerge output and to tophat output
cm_path=../cuffmerge_out/merged_asm/merged.gtf
th_path=../tophat_map

# create a directory for all the cuffquant results
mkdir ../data/cuffquant_out
cd ../data/cuffquant_out

# do the quantification
module load boost bowtie2 cufflinks

# iterate through all the samples
for d in ../tophat_map/*/
do
sn=$(basename "${d}")
echo "Quantifying" ${sn}
mkdir ${sample}
cd ${sample}
# change number of threads to suit your computer hardware
cuffquant -no-update-check -p 19 ${cm_path} ${th_path}/${sample}/accepted_hits.bam
cd ..
done > $log_file