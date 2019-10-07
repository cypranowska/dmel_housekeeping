#!/bin/bash
# Job name:
#SBATCH --job-name=hs_hisat2_map
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
#SBATCH --time=2-00:00:00
#
#SBATCH --mail-user=casegura@berkeley.edu
#SBATCH --mail-type=ALL
## Command(s) to run:

start=`date +%s`
echo "Loading hisat2"

module load hisat2
mkdir ../hisat2_map
cd ../hisat2_map

echo "Mapping reads with hisat2 in" $PWD

for f in ../fastq_trim/*1.trim.fastq.gz
do
base=$(basename "${f}")
sample=${base%_1.trim.fastq.gz}
echo "Aligning" ${sample}
hisat2 --dta -t -p 19 -x ../ref/dm6 -1 ../fastq_trim/${base} -2 ../fastq_trim/${sample}_2.trim.fastq.gz -S ${sample}.sam
done

end=`date +%s`
runtime=$((end-start))
echo "Done in" ${runtime}
