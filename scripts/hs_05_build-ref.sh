#!/bin/bash
# Job name:
#SBATCH --job-name=build_ref
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
#SBATCH --time=03:00:00
#
## Command(s) to run:

echo "Loading hisat2"

module load hisat2

cd ../ref

echo "Building hisat2 reference genome"

hisat2-build Drosophila_melanogaster.BDGP6.dna.toplevel.fa dm6
