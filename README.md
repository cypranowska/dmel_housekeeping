# dmel_housekeeping

## Description

This repository contains the bash and R scripts for downloading public modENCODE RNA-seq data from larval Drosophila melanogaster tissues, mapping the reads to the reference genome, quantifying the transcripts, and identifying genes that have well-correlated expression across multiple tissues at the same developmental stage. The overall goal of this project is to improve single-cell RNA-sequencing library normalization by identifying "housekeeping" genes that can serve as negative controls for generalized linear models to regress variation between samples while preserving real biological differences between cells.

This repository currently has scripts that download modENCODE data from NCBI SRA, performs quality control on the reads, aligns them to the Drosophila reference genome, normalizes the transcript abundances, and generates figures comparing gene expression in each sequenced tissue in the L3 larval data set. 

## Create output directories to add to .gitignore

In this project directory, make the appropriate output directories by running `mkdir -p data/fastq` and `mkdir R_out`. Add the `data` and `R_out` directories to .gitignore. Make sure the scripts are in the `scripts` directory.

## Scripts to be run (in order)

1. `hs_01_fastq_download.sh`
2. `hs_02_fastqc.sh`
3. `hs_03_trim.sh`
4. `hs_04_fastqc.sh`
5. `hs_05_build-ref.sh`
6. `hs_06_tophat.sh`
7. `hs_07_cufflinks.sh`
8. `hs_08_cuffmerge.sh`
9. `hs_09_cuffquant.sh`
10. `hs_10_cuffnorm.sh`

## Script details

### Fastq download
`hs_01_fastq_download.sh` uses SRA toolkit `fastq-dump` to fetch .fastq records pertaining to runs listed in `ref/srr_fastq_queue.csv`. Output is written to `data/fastq`.

### FastQC
`hs_02_fastqc.sh` and `hs_04_fastqc.sh` run `fastqc` on raw and trimmed reads respectively. Output is written to working directory and then moved to either `data/fastqc_raw` for the former and `data/fastqc_trimmed` for the latter. Fastq files are compressed with `gzip` in `hs_02_fastqc.sh` to reduce storage overhead.

### Trimmomatic
`hs_03_trim.sh` runs the Trimmomatic .jar file on each set of paired end reads with a sliding window of 4 and a quality threshold of 20. Script as it is written assumes at least 20 cores available on system. 

Output is written to `data/fastq` and trimmed reads (with `.trim.fastq.gz` extension) are moved to `data/fastq_trim`. 

### Building the reference genome


## System Requirements

This repository contains the entire computational pipeline for this project from data download from the NCBI Sequence Read Archive to visualizations of the normalized RNA-seq data in R. 

### Bash
Command-line tools and their versions used for this project are:

* SRA Toolkit ver 2.9.4 (CentOS Linux 64-bit distro)
* FastQC
* Trimmomatic ver 0.38
* Bowtie2 
* Samtools
* Tophat ver 2.1.1
* Cufflinks ver 2.2.1

### R
All R code was run on R version 3.5.2. The libraries required for R scripts are as follows:

* ggplot2
