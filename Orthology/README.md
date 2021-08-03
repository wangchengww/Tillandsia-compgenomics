
# Calling orthogroups between *T. fasciculata*, *T. leiboldiana* and *A.comosus*

Using orthofinder, I called orthogroups between three bromeliad gene model sets. The resulting orthogroups were later used to study synteny and gene family evolution.

# First run: all T. fasciulata and T. leiboldiana gene models included, regardless of location

The first orthofinder run involved all called gene models and was run to see the distribution of genes across the assembly, making a distinctions between one-to-one orthogroups and oorthogroups with other relationships. The idea behind this was that, if most gene models, especially gene models with one-to-one relationships are on main scaffolds, it confirms that the main scaffolds represent the chromosomes and we can limit downstream analyses to these sequences. Alternatively, if many genemodels are also present on smaller scaffolds, this may indicate fragmentation in our assembly and an anchoring step may be necessary.

**1) Input sequences**

The following peptide sequence datasets were provided for running orthofinder:

  - T. fasciculata: From our final annotation run (maker, work computer, fasciculata_EDTA_R1_full_prot_full_mrna_masked2.maker.output), after renaming all gene models in the protein file (see Annotation documentation). The file contains all gene models, including those that didn't blast in functional annotation: Tillandsia_fasciculata_v1.all.proteins.fasta (34886 sequences)
  - T.leiboldiana: From our final annotation run (maker, VSC4, leiboldiana_R1_alt.maker.output), after renaming all gene models. Again, the full set was used including sequences that didn't blast (38180)
  - Ananas comosus: From Ray Ming's website I downloaded the Acomosus F153 genome, which contains peptide sequences: http://www.life.illinois.edu/ming/LabWebPage/Downloads.html. I renamed the file to Acomosus_F153.20150427.proteins.fasta (27024)

**2) Running Orthofinder**

Orthofinder was run with version 2.4.1. with the following command:

    /apps/orthofinder/2.4.1/orthofinder \
	  -f /scratch/grootcrego/orthofinder/run_orthofinder_Tfas_Tlei_Acom \
	  -t 48

 **3) Compiling orthology and annotation results**

For the global run, I wanted to compile orthology information per gene with its location in the genome. First, I converted the N0.tsv file to a per-gene format which lists for each gene it's orthogroup and the orthogroup gene counts per species. Using the python script This was done by running the script `script_make_per_gene_og_table_calculate_counts.py` on N0.tsv.

NOTE: In previous runs I created the per-gene orthology table by combining the N0.tsv file with the Orthogroups.GeneCounts.tsv file. However, the counts in the latter file don't agree with the genes reported in the former since they stem from different approaches. So the newest version, counts are inmediately calculated from N0.tsv in the python script, since results from the Orthogroups/ folder are deprecated.

Since I am not interested in A.comosus specific orthogroups, these were removed from the new table:

    awk '!($4 == 0 && $5 == 0) {print $0}' \
    orthogroups_Tfas_Tlei_Acom.per_gene.txt > \
	orthogroups_Tfas-Tlei-Acom.no_Acom_specific_og.txt

This results in 20,507 orthogroups and 86,017 genes (31,551 Tfas genes, 32,568 Tlei genes, 21,898 Acom genes).

Then, I ran the script `script_compile_gff_info_og_table_all_species.py` to obtain the final per-gene table compiling location and orthology information:

	python2 compile_pergene_orthology_gff_ino.py \
	orthogroups_Tfas-Tlei-Acom.no_Acom_specific_og.txt \
	Tillandsia_leiboldiana_v1.2.edited_allfeatures.gff \
	Tillandsia_fasciculata_v1.2.edited_allfeatures.gff \
	Acom_annotation.txt orthogroups_Tfas-Tlei-Acom.full-set.no_Acom_specific_og.all_info.txt


**4) Analyzing the spatial and length distribution of global orthogroups**

With this final table, I investigated where orthologous gene are mostly found in the genome and how long they are. To do this, I separated the Tlei and Tfas genes in separate files:

    grep "Tlei" orthogroups_Tfas-Tlei-Acom.full-set.no_Acom_specific_og.all_info.txt > Tlei_orthology_info.txt
    grep "Tfas"  orthogroups_Tfas-Tlei-Acom.full-set.no_Acom_specific_og.all_info.txt > Tfas_orthology_info.txt

Lastly I added scaffolds lengths to these files and reshuffled the columns a bit to make it more logical for the next step with `add_scaffold_lengths.py`.

The actual analysis was done in the Rscript `Analyze_global_orthogroups.R`. The analysis showed that more than 99 % of 1-to-1 orthologs are on the main scaffolds in the case of both assemblies (> 1 Mb). This also led to the inclusion of the 26th scaffold in *T. leiboldiana*, as it contained a non-negligible amount of orthologs.

# Second run: Only gene models on main scaffolds (> 1 Mb) for T.fas and T. lei, all gene models of A.comosus

With the aim of studying synteny and gene family evolution between Tfas and Tlei, I decided to narrow down our list of genes to those which lie on the main scaffolds of the Tfas and Tlei assemblies. This way, we will avoid the presence of hidden TEs, as those are more likely to lie on rogue scaffolds.

**1) Input sequences and a bit of filtering**

To select the genes on those scaffolds, I made a list with the names of the main scaffolds `25_largest_scaffolds`:

    awk '$3 == "mRNA" {print $0}' Tillandsia_fasciculata_v1.2.edited.gff | grep -f 25_largest_scaffolds > \
	mRNA_entries_on_25_largest_scaffolds

Then, I selected these from the maker peptide file by extracting the feature ID and selecting the fasta-sequences matching these ID's:

    cut -f 9 mRNA_entries_on_25_largest_scaffolds | sed 's/;/\t/g' | cut -f 1 | sed 's/ID=//g' > \
    mRNA_entries_on_25_largest_scaffolds_IDonly
    seqkit grep -f mRNA_entries_on_25_largest_scaffolds_IDonly \
	Tillandsia_fasciculata_v1.all.proteins.fasta >  Tillandsia_fasciculata_v1.25chrom.proteins.fasta

I filtered out shorter isoforms (of which there were very few) and also peptide sequences < 40 amino acids.
The longest isoform was selected using the script `select_longest_isoform.py`.
Peptide sequences < 40 AA were filtered out with the following lines of code:

    awk '$2 > 40 {print $0}' \
    Tillandsia_fasciculata_v1.25chrom.longest_isoforms.proteins.fasta.fai \
	| cut -f 1 > proteins.w.min40AA
    seqkit grep -f proteins.w.min.40.AA \
	Tillandsia_fasciculata_v1.25chrom.longest_isoforms.proteins.fasta > \
	Tillandsia_fasciculata_v1.25chrom.longest_isoforms.min40AA.proteins.fasta

I moved these sequences to the orthofinder files in scratch under the new run folder run_orthofinder_Tfas_Tlei_Acom_25_scaffolds, where I renamed all protein files to easier names (T.fasciculata.fa, T.leiboldiana.fa, A.comosus.fa). Orthofinder was run with the same command as shown above.

**2) Compiling orthology with gff info and functional annotations, plus a bit more filtering**

This time, I did not remove the Acomosus information from the results file, though I did eliminate all orthogroups that have no orthologues in T.lei and T.fas (specific to A.comosus):
`awk '!($3 == 0 && $4 == 0) {print $0}' Orthogroups.GeneCount.tsv > orthogroup_counts_no_Acom_specific_og.txt`
This removed 553 orthogroups from the file. The orthogroup IDs in this file were then used to select non-pineapple specific orthogroups from N0.tsv.

IMPORTANT: because of the nature of the N0.tsv file, which contains nested HOG and does not agree with the files in the Orthogroups directory, I made significant changes to the script make_og_table.py. Previously, this table compiled information from N0.tsv (the genes in each orthogroup) with the information in Orthogroups.Genecounts.tsv (the counts of each orthogroup). After careful examination I realized these two files don't agree, since they stem from different approaches and the Orthogroups output is deprecated. Therefore, I generated the python script `script_make_per_gene_og_table_calculate_counts.mainscaffolds.py` that manually counted the number of genes per species in each orthogroup and appended this to the per-gene table

I then ran a more elaborate version of the compiling script which adds functional descriptions and gff file information for genes from all 3 species at once, named `script_compile_gff_info_with_og_tanle_all_species.mainscaffolds.py`:

    python2 ../make_big_og_table_all_info.py \
	orthogroups_no_Acom_specific_og.per_gene.txt Tillandsia_leiboldiana_v1.2.edited_allfeatures.gff \
	Tillandsia_fasciculata_v1.2.edited_allfeatures.gff Acom_annotation \
	orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.txt

I then searched for all orthogroups containing genes with TE description, by searching for the words:

    transposable
    transposon
    transposase
    Transposon
    Gag-Pol polyprotein
    Pro-Pol polyprotein
    virus

There were 1034 genes with these descriptions belonging to 264 orthogroups. I obtained the IDs of these orthogroups and then removed them from the table :

	grep -f tes orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.txt | cut -f 7 | \
	sort | uniq > orthogroups_containing_TES
    grep -v -f orthogroups_containing_TES \
	orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.txt > orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt

This removed a total of 6042 genes.
The file `orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt` is the final set of curated orthogroups that has been used in downstream analysis. It contains 70,963 genes and 19,101 orthogroups. Additional statistics on the orthology analysis can be found [here](https://docs.google.com/spreadsheets/d/1hv_Oe6MvV1fBuBkGIhvM1ihDwYviuY-LcYUGYKnIwNQ/edit#gid=324678736).
