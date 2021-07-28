# Studying synteny and rearrangements between *T. fasciculata* and *T. leiboldiana*

Using orthologous gene pairs, whole genome alignment and local alignments, I investigated synteny and identified rearrangements between both species and also *Ananas comosus*. This revealed a number of rearrangements and also provided insights into the genome architecture evolution in both species.

# Broad look at synteny with 1-to-1 orthologous pairs

I obtained all one-to-one orthologs from the file Orthologues/Orthologues_T.fasciculata/T.fasciculata__v__T.leiboldiana.tsv:
grep -v "," T.fasciculata__v__T.leiboldiana.tsv > one_to_one_orthologues. Using this file I can created the table for circlize as above, but with a modified script. This runs very slowly and should be optimized. Of the 14958 1-1 orthologues 75 were lost which were in orthogroups removed because they contained TEs. This resulted in the file circlize_table_one-to-one_orthology_Tfas-Tlei_25scaffolds.txt which I moved to my computer. There, I split the file into 25 scaffolds of Tfas with following bash script:

    i=1
    cat 25_largest_scaffolds | while read line ;do
     chrom=`echo "$line"`
     chrom_n=`echo "$line"| cut -f 2 -d'_'`
     echo $chrom
     echo $chrom_n
     echo $i
     awk -v pat="$chrom" '$2 ~ pat' circlize_table_one-to-one_orthology_Tfas-Tlei_25scaffolds.txt > chr${i}_${chrom_n}_Tfas_25chrom.txt;
     cut -f 1,2,3,4,6,7,8,9 chr${i}_${chrom_n}_Tfas_25chrom.txt > tmp
     mv tmp chr${i}_${chrom_n}_Tfas_25chrom.txt
     i=$((i+1))
    done

These files were used to rerun synteny plots.

I also curated the Gene Duplication results by extracting only the curated orthogroups in the Duplications.tsv files. This reduced the number of duplications significantly, especially in T. leiboldiana.

To study synteny, the above generated files of 1 to 1 orthologs per scaffolds were generated both per scaffold of Tfas and Tlei. Synteny plots were generated in both directions.
I also looked at which scaffolds within both species had most links of 1:1 orthologs, as synteny seems to be relatively well maintained. Generally, the majority of orthologs link a pair of chromosomes between Tfas and Tlei together. These pairs I now regard as syntenic and therefore "corresponding" chromosomes with probably common ancestry. I looked up these pairs by calculating the number of links per pairs:
`cut -f 2,6 chr1_2199_Tfas_25chrom.txt | sort | uniq -c`

Which gives:
    784 Scaffold_2199	Scaffold_10370
    3 Scaffold_2199	Scaffold_10423
    2 Scaffold_2199	Scaffold_10425
    2 Scaffold_2199	Scaffold_10429
    1 Scaffold_2199	Scaffold_10430
    1 Scaffold_2199	Scaffold_247

Here there is a clear indication that scaffold 2199 in Tfas (chr1) and scaffold 10370 (chr2) in Tlei correspond. The pairs of corresponding chromosomes between both species can be found in [this](https://docs.google.com/spreadsheets/d/1Gfj0WRwzEupbUZKON2OnO8psCsLD4j3tC6k3-Pz8hLs/edit#gid=0) spreadsheet on google drive.

# Whole genome alignment between Tfas, Tlei and A. comosus

Because this spreadsheet indicated a few potential rearrangements, I reran a whole genome alignment between Tfas and Tlei to obtain a second source of evidence:
    nucmer -p Tfas_Tlei_alignment_25chrom /proj/grootcrego/Genome_assemblies/fasciculata/4_final_assembly/tillandsia_fasciculata_assembly.sorted.25_scaffolds.fasta /proj/grootcrego/Genome_assemblies/leiboldiana/4_annotation/Tillandsia_leiboldiana_26_scaffolds.fasta

Both the synteny analysis using 1-1 orthologues and the dotplot indicated 3 possible rearrangements: 2 chromosome arm swaps and one merger of 2 chromosomes in T. leiboldiana.

To further assess this, I ran a nucmer and dotplot analysis between T.fasciculata / T. leiboldiana and A. comosus. These in fact showed the same rearrangements between Tlei and Acom and Tlei and Tfas.

# Local alignments between Tfas and Tlei using LastZ

I then performed LastZ alignments, which has a local blast strategy, to further confirm these breakpoits. For this, I reran hardmasking of both genomes with the newest version of TE annotation from EDTA. As discussed in "TE annotation", I only used the known TE families. The repeatmasker command was:

    RepeatMasker -dir . -gff -pa 4 -e ncbi -nolow -lib /proj/grootcrego/Genome_assemblies/fasciculata/4_final_assembly/TE_families_FULL_EDTA_Tfas_25chrom.fa /scratch/grootcrego/LastZ/Tillandsia_fasciculata_25_scaffolds.fasta

To run LastZ efficiently, I decided to always map the full T.lei genome to each T.fas chromosome and parallelize this. To create separate fasta files for each Tfas chromosome I used the following bash loop:

    cat tillandsia_fasciculata_chrnames.txt | while read line ; do  Name=`echo "$line"|awk '{print $1}'`;  echo $Name;  Replace=`echo "$line"|awk '{print $2}'`;  echo $Replace;  seqkit grep -p $Name Tillandsia_fasciculata_25_scaffolds.fasta.masked > Tfas_$Replace.fasta.masked;  sed  -i "s/$Name/Tfas_$Replace/g" Tfas_$Replace.fasta.masked; done

To replace chromosome names in the full T.lei genome, I ran:

    cat Tillandsia_leiboldiana_26_scaffolds_chrnames.txt | while read line ; do  Name=`echo "$line"|awk '{print $1}'`;  echo $Name;  Replace=`echo "$line"|awk '{print $2}'`;  echo $Replace;  sed  -i "s/$Name/Tlei_$Replace/g" Tillandsia_leiboldiana_26_scaffolds.fasta.masked; done

I ran LastZ alignments with a slurm array, see `run_lastz.sh`.

Then I used Tibo's scripts to filter alignments and visualize anchoring:

    for i in Tlei_vs_Tfas_chr*; do echo $i;  
     python2 convertMaftoCoordinates.py $i > ${i%.maf}.coord;
    done
    bash script_generate_matrix_alignblastz_results.sh test_out Tillandsia_leiboldiana_26_scaffolds.fasta.masked

The last command runs a bunch of scripts, among which a filtering step where the 95 % quantile for length and identity is determined, and anything under these thresholds is filtered out. The .filtered file was then used to visualize breakpoints using the rscript visualization_lastz_alignments.R

To remove noise, I wrote a python script that will eliminate all alignments which overlap with an alignment from a different chromosome above a certain threshold. This script was applied both on the length-filtered and non-filtered coord file:
    python2 script_filter_alignments_by_uniqness_threshold.py Tlei_vs_Tfas_allchrom_lastz.coord 0.9

 The LastZ alignments seem to support the rearrangements observed in dotplot and orthology analyses. In some cases, the breakpoints are clear-cut, in others not so much. I decided to zoom in on the clear breakpoints to find out a range in which the breakpoint may exist, and then look at the pacbio alignments in that region.
