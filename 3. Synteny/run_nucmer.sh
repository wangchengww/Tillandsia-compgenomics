#!/bin/bash
#
#SBATCH --job-name=Nucmer
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --partition=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<clara.groot.crego@univie.ac.at>

module load mummer

cd /scratch/grootcrego/nucmer

nucmer -p Tlei_Tfas_alignment_25chrom /scratch/grootcrego/LastZ/Tillandsia_leiboldiana_26_scaffolds.fasta /scratch/grootcrego/LastZ/Tillandsia_fasciculata_25_scaffolds.fasta
