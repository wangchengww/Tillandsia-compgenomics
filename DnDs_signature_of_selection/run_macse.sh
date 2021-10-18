#!/bin/bash
#
#SBATCH -J macse_ext
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_0384
#SBATCH --qos p71400_0384
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<clara.groot.crego@univie.ac.at>

list_of_orthologs="/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/dnds/filtering/orthogroups_final_subset_0.2_lengthdiff_1_complete.txt"
macse="/home/fs71400/grootcrego/software/macse_v2.05.jar"
fasta_dir="/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/dnds/fasta_files"
out_dir="/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/dnds/alignments/extensive_alignment_final_subset"

cat $list_of_orthologs | while read line; do
	orthogroup=`echo "$line"`
	java -jar $macse -prog alignSequences -seq $fasta_dir/${orthogroup}.fasta -local_realign_init 1 -local_realign_dec 1 -out_AA $outdir/${orthogroup}_AA.fasta -out_NT $outdir/${orthogroup}_NT.fasta;
done
