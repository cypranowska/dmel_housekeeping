#!/bin/bash
# this script normalizes the data using cuffnorm

log_file="$(date "+%m%d%Y-%H%M%S").log"

# create a variable with paths to cuffmerge output and to tophat output
cm_path=../cuffmerge_out/merged_asm/merged.gtf
cq_path=../cuffquant_out

# create a directory for first aggregated transcriptome
mkdir ../data/cuffnorm_out
cd ../data/cuffnorm_out

# create a variable array for each normalization option
norm[1]="classic-fpkm"
norm[2]="geometric"
norm[3]="quartile"

# do the merging
module load boost bowtie2 cufflinks
for i in 1 2 3
do
	norm_method=${norm[${i}]}
	echo "Normalizing with" $norm_method
	# change number of threads to suit computer hardware
	# cuffnorm input is formated so samples from the same tissue type (or biological condition) are separated by commas
	# and sample groups are separated by spaces
	cuffnorm -p 19 -library-norm-method ${norm_method} -o ${norm_method} ${cm_path} ${cq_path}/SRR070426/abundances.cxb,${cq_path}/SRR100269/abundances.cxb ${cq_path}/SRR070409/abundances.cxb,${cq_path}/SRR070410/abundances.cxb  ${cq_path}/SRR100268/abundances.cxb ${cq_path}/SRR070405/abundances.cxb,${cq_path}/SRR070406/abundances.cxb ${cq_path}/SRR070392/abundances.cxb,${cq_path}/SRR070393/abundances.cxb,${cq_path}/SRR111884/abundances.cxb,${cq_path}/SRR111885/abundances.cxb,${cq_path}/SRR350962/abundances.cxb,${cq_path}/SRR350963/abundances.cxb ${cq_path}/SRR070407/abundances.cxb,${cq_path}/SRR070408/abundances.cxb,${cq_path}/SRR070425/abundances.cxb
done > $log_file