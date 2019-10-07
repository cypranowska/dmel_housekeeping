#!/bin/bash
# Job name:
#SBATCH --job-name=featureCounts_hs
#
# Account:
#SBATCH --account=fc_nmjreg
#
# Partition:
#SBATCH --partition=savio
#
# Quality:
#
#SBATCH --qos=savio_normal
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=casegura@berkeley.edu
#
## Command(s) to run:

#creating an array to hold sample name

#create a variable with path to reference .gtf file
ref_path=/global/scratch/casegura/housekeeping/ref/Drosophila_melanogaster.BDGP6.95.gtf

#make directory for featureCounts output
cd ..
mkdir fc_out
cd fc_out

#load subread
module load subread

for f in ../hisat2_map/*.bam
do
sn=$(basename "${f}" .bam)
echo "Working on" $sn
mkdir ${sn}
cd ${sn}
featureCounts -p -T 19 -a ${ref_path} -o counts.txt ../../hisat2_map/${sn}.bam
cd ..
done
