#!/bin/bash
# this script merges the assembled transcripts of each sample into one merged transcriptome

log_file="$(date "+%m%d%Y-%H%M%S").log"

# create a variable with path to manifest .gtf file
manifest=../../ref/master-gtf.txt

# create a variable with path to reference .gtf file
ref_path=../../ref/Drosophila_melanogaster.BDGP6.95.gtf

# dreate a directory for the master aggregated transcriptome
mkdir ../data/cuffmerge_out
cd ../data/cuffmerge_out

# do the merging
module load boost bowtie2 cufflinks
cuffmerge -g $ref_path $manifest > $log_file