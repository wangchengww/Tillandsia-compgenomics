#!/bin/bash
#
#SBATCH -J macse_121
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_0096
#SBATCH --qos p71400_0096
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<clara.groot.crego@univie.ac.at>

list_of_orthologs="/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/dnds/dnds_paralogs/orthologs-121.perOG.txt"
macse="/home/fs71400/grootcrego/software/macse_v2.05.jar"
fasta_dir="/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/dnds/dnds_paralogs/fasta_seq_per_OG_121"

cd alignments_121/
cat $list_of_orthologs | while read line; do
	orthogroup=`echo "$line"|awk '{print $1}'`
	java -jar $macse -prog alignSequences -seq $fasta_dir/${orthogroup}_Tfas-copy1.fasta -local_realign_init 1 -local_realign_dec 1 ;
	java -jar $macse -prog alignSequences -seq $fasta_dir/${orthogroup}_Tfas-copy2.fasta -local_realign_init 1 -local_realign_dec 1 ;
done
