#!/bin/bash
#
#SBATCH -J featureCounts
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_0384
#SBATCH --qos p71400_0384
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<clara.groot.crego@univie.ac.at>

gff=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/Tfas_assembly/Tillandsia_fasciculata_v1.2.edited_allfeatures.25chrom.gff
output=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/counts.Tfas_Tlei_6_timepoints.txt
reference=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/Tfas_assembly/Tillandsia_fasciculata_25_scaffolds.fasta
input=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/2_mapped/*.bam

mkdir tmp

featureCounts -a $gff -o $output -g ID -t gene -G $reference -T 48 -p -s 2 --tmpDir tmp/ $input
