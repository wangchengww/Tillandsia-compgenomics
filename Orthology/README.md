
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

For the global run, I wanted to compile orthology information per gene with its location in the genome. This was done first by


 The main information I was interested in from the orthofinder run, was to which orthogroup each gene belonged, what the gene count for each species was in each orthogroup, and what the functional annotations of each orthogroup are. I compiled all this information into one per-gene table where a line represents one gene and contains its coordinates, its orthogroup, the orthogroup counts for each species, and its functional annotation. This table contains genes from all three species.

I obtained this table with the python scripts


  To investigate all of this, I decided to create a table for each species parting from orthofinder's results and the assembly's gff file. The table will contain a gene model on each line, with scaffold, start and end position, length, orthogroup and number of genes of each species in the orthogroup. This way we will be able to compare one-to-one orthologues to one-to-many orthologues parting from the same table.

  First, I made a table containing orthogroup name and number of genes of each species parting from the orthogroup output file Orthogroups.GeneCount.tsv in /scratch/grootcrego/orthofinder/run_orthofinder_Tfas_Tlei_Acom/OrthoFinder/Results_Nov02/Orthogroups on the cube
  NOTE: after the running orthofinder 2.4.0, the file used for this computation was: /scratch/grootcrego/orthofinder/run_orthofinder_Tfas_Tlei_Acom/OrthoFinder/Results_Nov18/Orthogroups/Orthogroups.GeneCount.tsv
    cut -f 1,3,4 Orthogroups.GeneCount.tsv | awk '!($2 == 0 && $3 == 0) {print $0}' > orthogroup_counts_Tfas_Tlei.txt
  This only selects information on Tfas and Tlei and removes all orthogroups that were unique to Acom (0 genes in Tfas and Tlei). This resulted in 20,507 orthogroups with at least one gene for one of the two species per group.
  I then extracted the names of these orthogroups, removing the first line:
    cut -f 1 orthogroup_counts_Tfas_Tlei.txt | tail -n+2 > tmp
  With this list, I extracted the names of genes per relevant orthogroup from Orthogroups.txt
  NOTE: after the running orthofinder 2.4.0, the file used for this computation was: /scratch/grootcrego/orthofinder/run_orthofinder_Tfas_Tlei_Acom/OrthoFinder/Results_Nov18/Phylogenetic_Hierarchical_Orthogroups/N0.tsv
    grep -w -f tmp Orthogroups.txt > orthogroups_Tfas_Tlei.txt

  I then converted the orthogroups file into a table where each line was a gene and its corresponding orthogroup. This I did with the python script make_og_per_gene.py (edited for run2):
    #!/usr/bin/env python
    import sys
    ogroups = open(sys.argv[1])
    outputfilename=sys.argv[1].replace(".txt",".per_gene.txt")
    output=open(outputfilename,'w')

    for line in ogroups:
       line = line.replace('\r\n','')
       splitted_line = line.split('\t')
       OG_ID = splitted_line[1]
       genes_in_OG = splitted_line[4]+", "+splitted_line[5]
       genes_in_OG2 = genes_in_OG.split(", ")
       for gene in genes_in_OG2:
           if gene == '':
               continue
           else:
               line_to_print = gene+"\t"+OG_ID+"\n"
               output.write(line_to_print)

Then I merged this table into the counts table, to obtain one final table where each row is a gene with its corresponding orthogroup and the count of Tfas and Tlei genes in this respective orthogroup. This was done with the python script make_og_table.py:

    #!/usr/bin/env python

    import sys

    per_gene = open(sys.argv[1])
    counts = open(sys.argv[2])
    outputfilename=sys.argv[1].replace(".txt",".table.txt")
    output=open(outputfilename,'w')

    og_dict = {}  # I initiate the dictionary
    for line1 in counts.readlines()[1:]:
	   line1 = line1.replace('\n','') # rm the return carriage
	   splitted_line1 = line1.split('\t', 1) # split the line only at first occurring tab
	   og = splitted_line1[0]
	   counts_line1 = splitted_line1[1]
     og_dict[og] = counts_line1

    for line2 in per_gene.readlines():
     line2 = line2.replace('\n','') # rm the return carriage
     splitted_line2 = line2.split('\t')
     og_id = splitted_line2[1]
     if (og_dict.has_key(og_id)):
        line_to_print = line2+"\t"+og_dict[og_id]+"\n"
        output.write(line_to_print)

The resulting table looks like this:
  Tfasc_v1.31994-RA	OG0000000	1	2100
  Tlei_v1.00256-RA	OG0000000	1	2100
  Tlei_v1.00257-RA	OG0000000	1	2100
  Tlei_v1.00262-RA	OG0000000	1	2100
  Tlei_v1.00273-RA	OG0000000	1	2100
  Tlei_v1.00325-RA	OG0000000	1	2100
  Tlei_v1.00337-RA	OG0000000	1	2100
  Tlei_v1.00349-RA	OG0000000	1	2100
  Tlei_v1.00362-RA	OG0000000	1	2100
  Tlei_v1.00363-RA	OG0000000	1	2100

I then ran the python script make_big_og_table.py to combine information from each assembly's gff and the above table:

  #!/usr/bin/env python

  import sys

  orthogroups = open(sys.argv[1])
  gff = open(sys.argv[2])
  outputfilename=sys.argv[3]
  output=open(outputfilename,'w')

  # 1st, make a dictionnary of the gff with the information we want to transfer to the final table (scaffold, end and start position)

  gff_dict = {}
  for line1 in gff.readlines():
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t') # split the line regarding the tabulations
    if splitted_line1[2] == "mRNA":
        pre_ID = splitted_line1[8]
        pre_ID = pre_ID.split(";")
        ID = pre_ID[0]
        ID = ID.replace('ID=', '')
        info = splitted_line1[0]+"\t"+splitted_line1[3]+"\t"+splitted_line1[4]
        gff_dict[ID]=info
    else:
        continue

  # 2nd, iterate over the orthogroup list and print each line together with the gff information into a new table

  for line2 in orthogroups.readlines():
    line2 = line2.replace('\n','') # rm the return carriage
    splitted_line2 = line2.split('\t') # split the line regarding the tabulations
    ID_OG = splitted_line2[0]
    if gff_dict.has_key(ID_OG):
        line_to_print = ID_OG+"\t"+gff_dict[ID_OG]+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\n"
        output.write(line_to_print)
    else:
        line_to_print = ID_OG+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\n"
        output.write(line_to_print)

 This was run with the following commands:
  python2 make_big_og_table.py orthogroups_Tfas_Tlei.per_gene.table.txt /proj/grootcrego/Genome_assemblies/fasciculata/4_final_assembly/Tillandsia_fasciculata_v1.2.edited_allfeatures.gff Tfas_orthology_info_per_scaffold.txt
  python2 make_big_og_table.py orthogroups_Tfas_Tlei.per_gene.table.txt /proj/grootcrego/Genome_assemblies/leiboldiana/4_annotation/Tillandsia_leiboldiana_v1.2.edited_allfeatures.gff Tlei_orthology_info_per_scaffold.txt

 I then removed all Tfas genes from the Tlei file and viceversa with grep:
  grep -v "Tlei" Tfas_orthology_info_per_scaffold.txt > tmp
  mv tmp Tfas_orthology_info_per_scaffold.txt

The table looks like this:
Tlei_v1.00256-RA	Scaffold_8399	3535789	3536604	OG0000000	1	2100
Tlei_v1.00257-RA	Scaffold_8399	3544423	3545327	OG0000000	1	2100
Tlei_v1.00262-RA	Scaffold_8399	3653612	3664165	OG0000000	1	2100
Tlei_v1.00273-RA	Scaffold_8399	3959136	3967127	OG0000000	1	2100
Tlei_v1.00325-RA	Scaffold_8399	5112048	5118927	OG0000000	1	2100
Tlei_v1.00337-RA	Scaffold_8399	5899662	5939490	OG0000000	1	2100
Tlei_v1.00349-RA	Scaffold_8399	6636290	6636577	OG0000000	1	2100
Tlei_v1.00362-RA	Scaffold_8399	7359129	7375304	OG0000000	1	2100
Tlei_v1.00363-RA	Scaffold_8399	7419222	7419791	OG0000000	1	2100
Tlei_v1.00364-RA	Scaffold_8399	7476401	7477620	OG0000000	1	2100

With this table, we can investigate the number of orthologous genes per scaffold and the proportion of 1:1 to 1:many to have an idea of gene distribution. I also added scaffold lengths to each row with the script add_scaffold_lengths_to_og_table.py. These tables were loaded in R for that analysis.

Additionally, I created a table containing the one-to-one orthologues per line along with their corresponding scaffolds and positions. This file was used to create the circular plot with circlize in R. To generate this file, I first selected just one_to_one orthologues:
  awk '$7 == 1 && $8 == 1 {print $0}' Tfas_orthology_info_per_scaffold.run2.txt > Tfas_one-to-one_orthologues.txt
  awk '$7 == 1 && $8 == 1 {print $0}' Tlei_orthology_info_per_scaffold.run2.txt > Tlei_one-to-one_orthologues.txt

Then I created the table (circlize_table_one-to-one_orthology_Tfas-Tlei.txt) with a homemade python script:
  python2 make_orthology_table_for_circlize.py Tfas_one-to-one_orthologues.txt Tlei_one-to-one_orthologues.txt
This table was transferred to my PC and split into each chromosome of T.fasciculata. The script for building the links in the circular plot is called circlize.R


2.3. Second run, genes only on main scaffolds
---

With the aim of studying gene family evolution between Tfas and Tlei, we decided to narrow down our list of genes to those which lie on the 25 / 26 largest chromosomes. This way, we will avoid the presence of genes which may be duplicated or hidden TEs, as those are more likely to lie on rogue scaffolds.
To select the genes on those scaffolds:
  awk '$3 == "mRNA" {print $0}' Tillandsia_fasciculata_v1.2.edited.gff | grep -f 25_largest_scaffolds > mRNA_entries_on_25_largest_scaffolds
Then, I selected these from the protein file:
  cut -f 9 mRNA_entries_on_25_largest_scaffolds | sed 's/;/\t/g' | cut -f 1 | sed 's/ID=//g' > tmp
  mRNA_entries_on_25_largest_scaffolds_IDonly
  seqkit grep -f mRNA_entries_on_25_largest_scaffolds_IDonly Tillandsia_fasciculata_v1.all.proteins.fasta > Tillandsia_fasciculata_v1.25chrom.proteins.fasta

Then, I selected the longest isoform of each gene using the script select_longest_isoform.py. The resulting file was then indexed with samtools to look at sizes. I decided to remove all protein sequences with less than 40 amino acids (642):
awk '$2 > 40 {print $0}' Tillandsia_fasciculata_v1.25chrom.longest_isoforms.proteins.fasta.fai | cut -f 1 > proteins.w.min40AA
seqkit grep -f proteins.w.min.40.AA Tillandsia_fasciculata_v1.25chrom.longest_isoforms.proteins.fasta > Tillandsia_fasciculata_v1.25chrom.longest_isoforms.min40AA.proteins.fasta

I moved these sequences to the orthofinder files in scratch under the new run folder run_orthofinder_Tfas_Tlei_Acom_25_scaffolds, where I renamed all protein files to easier names (T.fasciculata.fa, T.leiboldiana.fa, A.comosus.fa).

The results of this run are in: /scratch/grootcrego/orthofinder/run_orthofinder_Tfas_Tlei_Acom_25_scaffolds/OrthoFinder/Results_Jan15_1/

This time, I did not remove the Acomosus information from the results file, though I did eliminate all orthogroups that have no orthologues in T.lei and T.fas (specific to A.comosus):
  awk '!($3 == 0 && $4 == 0) {print $0}' Orthogroups.GeneCount.tsv > orthogroup_counts_no_Acom_specific_og.txt
This removed 553 orthogroups from the file.

IMPORTANT: because of the nature of the N0.tsv file, which contains nested HOG and does not agree with the files in the Orthogroups directory, I made significant changes to the script make_og_table.py. Previously, this table complied information from N0.tsv (the genes in each orthogroup) with the information in Orthogroups.Genecounts.tsv (the counts of each orthogroup). After careful examination I realized these two files don't agree, since they stem from different approaches and the Orthogroups output is deprecated. Therefore, I generated a python script that manually counted the number of genes per species in each orthogroup and appended this to the per.gene table:

import sys
ogroups = open(sys.argv[1])
outputfilename=sys.argv[1].replace(".txt",".per_gene.txt")
output=open(outputfilename,'w')
og_dict = {}
for line in ogroups:
    line = line.replace('\r\n','')
    splitted_line = line.split('\t')
    OG_ID = splitted_line[1]
    genes_in_OG =splitted_line[3]+", "+splitted_line[4]+", "+splitted_line[5]
    genes_in_OG2 = genes_in_OG.split(", ")
    if OG_ID not in og_dict.keys():
        og_dict[OG_ID] = genes_in_OG
    else:
        og_dict[OG_ID] = og_dict[OG_ID]+", "+genes_in_OG

for key in og_dict:
    OG_ID = key
    genes_in_OG2 = og_dict[key]
    genes_in_OG2 = genes_in_OG2.split(", ")
    Acom = 0
    Tfas = 0
    Tlei = 0
    for gene in genes_in_OG2:
        if gene.startswith("Aco"):
            Acom = Acom + 1
        if gene.startswith("Tfasc"):
            Tfas = Tfas + 1
        if gene.startswith("Tlei"):
            Tlei = Tlei + 1
    for gene in genes_in_OG2:
        if gene == '':
            continue
        else:
            Acom = str(Acom)
            Tfas = str(Tfas)
            Tlei = str(Tlei)
            line_to_print = gene+"\t"+OG_ID+"\t"+Acom+"\t"+Tfas+"\t"+Tlei+"\n"
            print(line_to_print)

 This produces a table with the correct number of genes (rows), 77,005.
 After doing this the python script make_og_table.py becomes obsolete. I then ran a more elaborate version of make_big_og_table.py which adds functional descriptions and gff file information for genes from all 3 species:

  #!/usr/bin/env python
  import sys
  orthogroups = open(sys.argv[1])
  gff_Tlei = open(sys.argv[2])
  gff_Tfas = open(sys.argv[3])
  gff_Acom = open(sys.argv[4])
  outputfilename=sys.argv[5]
  output=open(outputfilename,'w')

  # 1st, make a dictionnary of the gff with the information we want to transfer to the final table (scaffold, end and start position)

  gff_dict_Tlei = {}
  for line1 in gff_Tlei.readlines():
     line1 = line1.replace('\n','') # rm the return carriage
     splitted_line1 = line1.split('\t') # split the line regarding the tabulations
     if splitted_line1[2] == "mRNA":
         pre_ID = splitted_line1[8]
         pre_ID = pre_ID.split(";")
         ID = pre_ID[0]
         ID = ID.replace('ID=', '')
         description_line = pre_ID[7]
         GO = "NA"
         for element in pre_ID:
             if element.startswith("Ontology_id"):
                 ont = element
                 ont_split = ont.split('=')
                 GO = ont_split[1]
                 info = splitted_line1[0]+"\t"+splitted_line1[3]+"\t"+splitted_line1[4]+"\t"+GO+"\t"+description_line
                 gff_dict_Tlei[ID]=info
             else:
                 continue

 gff_dict_Tfas = {}
for line1 in gff_Tfas.readlines():
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t') # split the line regarding the tabulations
    if splitted_line1[2] == "mRNA":
        pre_ID = splitted_line1[8]
        pre_ID = pre_ID.split(";")
        ID = pre_ID[0]
        ID = ID.replace('ID=', '')
        description_line = pre_ID[7]
        GO = "NA"
        for element in pre_ID:
            if element.startswith("Ontology_id"):
                ont = element
                ont_split = ont.split('=')
                GO = ont_split[1]
        info = splitted_line1[0]+"\t"+splitted_line1[3]+"\t"+splitted_line1[4]+"\t"+GO+"\t"+description_line
        gff_dict_Tfas[ID]=info
    else:
        continue

gff_dict_Acom = {}
for line1 in gff_Acom.readlines():
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t') # split the line regarding the tabulations
    ID = splitted_line1[0]
    description_line = splitted_line1[3]
    GO = "NA"
    for element in splitted_line1:
        if element.startswith("GO"):
            GO = element
    position_line = splitted_line1[1]
    position_line = position_line.replace(':', '-')
    position_line = position_line.split('-')
    info = position_line[0]+"\t"+position_line[1]+"\t"+position_line[2]+"\t"+GO+"\t"+description_line
    gff_dict_Acom[ID]=info


# 2nd, iterate over the orthogroup list and print each line together with the gff information into a new table

for line2 in orthogroups.readlines():
    line2 = line2.replace('\n','') # rm the return carriage
    splitted_line2 = line2.split('\t') # split the line regarding the tabulations
    ID_OG = splitted_line2[0]
    if gff_dict_Tlei.has_key(ID_OG):
        line_to_print = ID_OG+"\t"+gff_dict_Tlei[ID_OG]+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\t"+splitted_line2[4]+"\n"
        output.write(line_to_print)
    elif gff_dict_Tfas.has_key(ID_OG):
        line_to_print = ID_OG+"\t"+gff_dict_Tfas[ID_OG]+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\t"+splitted_line2[4]+"\n"
        output.write(line_to_print)
    elif gff_dict_Acom.has_key(ID_OG):
        line_to_print = ID_OG+"\t"+gff_dict_Acom[ID_OG]+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\t"+splitted_line2[4]+"\n"
        output.write(line_to_print)

 This script was named make_big_og_table_all_info.py and run in the following command:
   python2 ../make_big_og_table_all_info.py orthogroups_no_Acom_specific_og.per_gene.txt /proj/grootcrego/Genome_assemblies/leiboldiana/4_annotation/Tillandsia_leiboldiana_v1.2.edited_allfeatures.gff /proj/grootcrego/Genome_assemblies/fasciculata/4_final_assembly/Tillandsia_fasciculata_v1.2.edited_allfeatures.gff /proj/grootcrego/genomic_resources/Acomosus_resource/tmp orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.txt

I then searched for all orthogroups containing genes with TE description, by searching for the words:
transposable
transposon
transposase
Transposon
Gag-Pol polyprotein
Pro-Pol polyprotein
virus
There were 1034 genes with these descriptions belonging to 264 orthogroups. I obtained the IDs of these OG and then removed them from the table :
   grep -f tes orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.txt | cut -f 7 | sort | uniq > orthogroups_containing_TES
   grep -v -f orthogroups_containing_TES orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.txt > orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt
This removed a total of 5687 genes.

With this curated file, I calculated the stats in the google drive file.

One important distinction that Orthofinder makes is between orthogroups and orthologs. An orthogroup is a group of genes descended from a single gene in the LCA of that group of species. It will include all gene duplications that occurred along the way. Orthologs are pairs of genes descended of a single gene in a common ancestor. So, if we have an orthogroup of 2 genes in Tlei and 3 genes in Tfas, Orthofinder can estimate which gene duplicated. One gene in Tlei may then have a one-to-one relationship to another gene in Tfas, and the other Tlei gene may have a one-to-two relationship to the remaining two Tfas genes.
The number of orthogroups with just 1 gene in T.lei and 1 gene in T.fas is about 13,128. Yet Orthofinder discovered about 14,963 one-to-one orthologues. These additional 1800 one-to-one genes actually belong to orthogroups that contain paralogs in one or the other species (in other words, at some point between the LCA of Tlei and Tfas and now, gene duplication occurred).
The ortholog relationships are useful for studies of gene duplication and family expansion. Therefore, in those analyses I would work on the ortholog level. For analyses of synteny, the larger group (orthologue level) should also suffice, so I decided to rerun that with the 14000 genes.
