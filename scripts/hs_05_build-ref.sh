#!/bin/bash
# this script is for building and indexing the Dmel 6 reference genome prior to mapping reads

log_file="$(date "+%m%d%Y-%H%M%S").log"

echo "Loading bowtie2"
module load bowtie2 samtools

cd ../ref

echo "Building bowtie2 reference genome"
bowtie2-build Drosophila_melanogaster.BDGP6.dna.toplevel.fa dm6 > $log_file

echo "Indexing reference fasta file"
samtools faidx Drosophila_melanogaster.BDGP6.dna.toplevel.fa >> $log_file