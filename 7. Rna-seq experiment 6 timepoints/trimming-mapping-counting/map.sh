#!/bin/bash
#
#SBATCH -J STAR
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_0096
#SBATCH --qos p71400_0096
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<clara.groot.crego@univie.ac.at>


wd=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/1_trimmed/
ref_dir=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/Tfas_assembly
ref=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/Tfas_assembly/Tillandsia_fasciculata_25_scaffolds.fasta
gff_file=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/Tfas_assembly/Tillandsia_fasciculata_v1.2.edited_allfeatures.25chrom.gff
star=/home/fs71400/grootcrego/software/STAR-2.7.9a/source/STAR

cd $wd

# Index Genome
#$star --runMode genomeGenerate --genomeDir $ref_dir --genomeFastaFiles $ref --sjdbGTFfile $gff_file --sjdbGTFtagExonParentTranscript Parent  --sjdbOverhang 149 --runThreadN 24 --genomeSAindexNbases 13

R1=($wd/Tlei*.pair1.truncated)
R2=($wd/Tlei*.pair2.truncated)
for ((i=0;i<=${#R1[@]};i++))
do
 $star --genomeDir $ref_dir --readFilesIn "${R1[i]}" "${R2[i]}" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "${R1[i].%.pair1.truncated}_toTfas." --limitBAMsortRAM 10240000000 --runThreadN 48
 mv *.bam *Log* /gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/2_mapped/
done
