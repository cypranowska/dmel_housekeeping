#!/bin/bash
# this script downloads SRR fastqs listed in srr_fastq_queue.csv and puts them data/fastq 

log_file="$(date "+%m%d%Y-%H%M%S").log"
outdir="../data/fastq"
infile="../ref/srr_fastq_queue.csv"

for SRRID in $(tail -n+7 $infile | cut -d"," -f1); do

echo $SRRID | tee -a "$log_file"
if [ ! -e $outdir'/'$SRRID'_1.fastq' ] | [ ! -e $outdir'/'$SRRID'_2.fastq' ]
then
~/modules/sratoolkit.2.9.4-1-centos_linux64/bin/fastq-dump --outdir $outdir --split-files -I $SRRID | tee -a "$log_file" && rm $HOME/ncbi/public/sra/*.sra* && echo "successfully downloaded "$SSRID | tee -a "$log_file"
fi
done
