#!/bin/bash
#
#SBATCH --job-name=lastz
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH -a 1-25%5
#SBATCH --partition=basic
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<clara.groot.crego@univie.ac.at>

module load lastz/1.04.03
lastz Tlei_chr${SLURM_ARRAY_TASK_ID}.fasta.masked Tillandsia_fasciculata_25_scaffolds.fasta.masked --notransition --step=10 --gapped --chain --gfextend --format=maf > Tfas_vs_Tlei_chr${SLURM_ARRAY_TASK_ID}_LASTZ_alignment.maf
