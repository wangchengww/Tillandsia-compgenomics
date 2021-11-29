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


wd=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/1_trimmed/subsample_mappingtest
ref_dir=/home/fs71400/grootcrego/Acomosus_resources/Acomosus_reference_genome
ref=/home/fs71400/grootcrego/Acomosus_resources/Acomosus_reference_genome/Acomosus_ASM154086v1_2016.fna
gff_file=/home/fs71400/grootcrego/Acomosus_resources/Acomosus_reference_genome/GCF_001540865.1_ASM154086v1_genomic.gff.gz
star=/home/fs71400/grootcrego/software/STAR-2.7.9a/source/STAR

cd $wd

# Index Genome
$star --runMode genomeGenerate --genomeDir $ref_dir --genomeFastaFiles $ref --sjdbGTFfile $gff_file --sjdbGTFtagExonParentTranscript Parent  --sjdbOverhang 149 --runThreadN 24 --genomeSAindexNbases 13

R1=($wd/*.pair1.truncated.subsample)
R2=($wd/*.pair2.truncated.subsample)
for ((i=0;i<=${#R1[@]};i++))
do
 $star --genomeDir $ref_dir --readFilesIn "${R1[i]}" "${R2[i]}" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "${R1[i].%.pair1.truncated.subsample}_toAco" --limitBAMsortRAM 10240000000 --runThreadN 48
 mv *.bam *Log* /gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/2_mapped/subsample_mappingtest
done
