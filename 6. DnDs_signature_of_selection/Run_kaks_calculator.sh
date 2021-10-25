#!/bin/bash
#
#SBATCH -J KaKs_Calculator
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks-per-core=1
#SBATCH -p mem_0384
#SBATCH --qos p71400_0384
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<clara.groot.crego@univie.ac.at>

modes=("stringent_filtering" "relaxed_reldiff" "relaxed_comp" "relaxed_reldiff_comp")
convert_script="/home/fs71400/grootcrego/software/AlignmentProcessor/bin/03_ConvertFasta.py"
kaks_script="/home/fs71400/grootcrego/software/AlignmentProcessor/bin/04_CallKaKs.py"
input="/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/dnds/alignments"
output="/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/dnds/dnds_calculations/KaKs_calculator"

for i in ${modes[*]}; do
	echo $i
	python $convert_script -i $input/extensive_alignments_$i -o $output/extensive_$i --axt
	python $kaks_script -i $output/extensive_$i -o $output/extensive_$i -m MYN
	mv $output/extensive_$i/KaKs.csv $output/extensive_$i/KaKs_extensive_$i.txt
done
