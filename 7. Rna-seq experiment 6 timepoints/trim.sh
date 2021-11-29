#!/bin/bash
#
#SBATCH -J Tlei_Qtrim
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_0384
#SBATCH --qos p71400_0384
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<clara.groot.crego@univie.ac.at>

input_dir='/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/0_raw'
output_dir='/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/1_trimmed'

R1=($input_dir/Tlei*.R1.fastq.gz)
R2=($input_dir/Tlei*.R2.fastq.gz)
for ((i=0;i<=${#R1[@]};i++))
do
  AdapterRemoval --file1 "${R1[i]}" --file2 "${R2[i]}" --basename "${R1[i]%.R1.fastq.gz}" --trimns --trimqualities --minquality 20 --trimwindows 12 --minlength 36
  fastqc "${R1[i]%.R1.fastq.gz}.pair1.truncated" "${R2[i]%.R2.fastq.gz}.pair2.truncated"
  mv *truncated *discarded *settings $output_dir
  mv *fastqc* $output_dir/qc
done
