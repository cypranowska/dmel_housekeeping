#!/bin/bash
# Job name:
#SBATCH --job-name=hs_sort_map
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
echo "Loading samtools"

module load samtools
cd ../hisat2_map

echo "Sorting .sam files in" $PWD

for f in *.sam
do
base="${f%.*}"
echo "Sorting" ${f}
samtools sort -@ 19 -o ${base}.bam ${f}
done

echo "Removing old .sam files"
rm *.sam

end=`date +%s`
runtime=$((end-start))
echo "Done in" ${runtime}
