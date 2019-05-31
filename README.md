# dmel_housekeeping

## Description

This repository contains the bash and R scripts for downloading public modENCODE RNA-seq data from larval Drosophila melanogaster tissues, mapping the reads to the reference genome, quantifying the transcripts, and identifying genes that have well-correlated expression across multiple tissues at the same developmental stage. The overall goal of this project is to improve single-cell RNA-sequencing library normalization by identifying "housekeeping" genes that can serve as negative controls for generalized linear models to regress variation between samples while preserving real biological differences between cells.

This repository currently has scripts that download modENCODE data from NCBI SRA, performs quality control on the reads, aligns them to the Drosophila reference genome, normalizes the transcript abundances, and generates figures comparing gene expression in each sequenced tissue in the L3 larval data set. 

## Create output directories to add to .gitignore

In this project directory, make the appropriate output directories by running `mkdir -p bash_out/fastq_raw` and `mkdir R_out`. Add the `data`, `bash_out`, and `R_out` directories to .gitignore. Make sure the scripts are in the scripts directory. 

