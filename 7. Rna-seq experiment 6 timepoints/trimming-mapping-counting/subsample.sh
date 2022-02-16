#!/bin/bash
#
#SBATCH -J subsample
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_0384
#SBATCH --qos p71400_0384
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<clara.groot.crego@univie.ac.at>

input_dir='/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/1_trimmed/subsample_mappingtest'


R1=($input_dir/*.pair1.truncated)
R2=($input_dir/*.pair2.truncated)
for ((i=0;i<=${#R1[@]};i++))
do
  seqkit sample -p 0.1 -s 11 -o "${R1[i]}.subsample" ${R1[i]}
  seqkit sample -p 0.1 -s 11 -o "${R2[i]}.subsample" ${R2[i]}
done
