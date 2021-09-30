# Inferring signatures of selection by calculating Dn/Ds ratios between T. fasciulata and T. leiboldiana

In this section, we infer Dn / Ds ratios for all single-copy pairs of orthologous genes between T. fasciculata and T. leiboldiana. Part of this work was done by Francesca Beclin. We obtained Dn / Ds ratios by aligning the orthologous sequences using a codon-aware aligner (MACSE) and inferring ratios with KaKs_Calculator and PAML.

# Selecting single-copy orthogroups and compiling per-orthogroup fastafiles

First I compiled the sequences of each gene into per-orthogroup fastafiles. For this, I extracted the one-to-one genes obtained in orthofinder analyses from the annotation gff files of each assembly. The orthofinder file used for this was: /scratch/grootcrego/orthofinder/run_orthofinder_Tfas_Tlei_Acom_25_scaffolds/OrthoFinder/Results_Jan22/Phylogenetic_Hierarchical_Orthogroups/one-to-one_orthogroups_Tfas_Tlei.txt and was created in the postprocessing of orthofinder results (see Orthology&Synteny notes). There are 13,128 orthologous pairs in this file.

To extract the features belonging only to one-one orthologous genes:

    grep -f one-to-one_orthogroups_Tfas_Tlei.IDonly.txt \
    Tillandsia_fasciculata_v1.2.edited_allfeatures.gff > \
    Tillandsia_fasciculata_v1.2.edited_allfeatures.one-to-one_orthologs.gff

Then I extracted the coding sequence features and replaced ID with Name to make the below python script work:

    awk '$3 == "CDS" {print $0}' [gff] > [gff.CDSonly]
    sed 's/ID=/Name=/g' [gff.CDSonly] > [gff.CDSonly2]

I then modified the scaffold names in the gff file from the old system ("scaffold_8339") to a new system ("chr1"):

    cat chrnames_Tfas.txt | while read line ;
    do
     Name=`echo "$line"|awk '{print $1}'`;
     echo $Name;
     Replace=`echo "$line"|awk '{print $2}'`;
     echo $Replace;
     sed  -i "s/$Name/$Replace/g" Tillandsia_fasciculata_v1.2.one-to-one_orthologs.CDSonly.gff;
    done

Then I extracted the fasta sequences for each orthologue and species using with `script_CutSeq_modified.py`, which is a modified form of one of Tibo's scripts to include a more informative header:

    while read line; do python2 cutSeqGff_mod.py  \
    Tfas_assembly/per_chrom_fasta/Tfas_$line.fasta \
    Tillandsia_fasciculata_v1.2.one-to-one_orthologs.CDSonly.gff \
    $line CDS; done < chr_Tfas.txt

I also compiled a list of orthologous gene pairs per orthogroup, which is used to compile the fasta sequences. This was done with the python script `script_compile_per-orthogroup_gene_list.py`

Then I used the following bash script and the resulting file to compile orthologous fasta sequences into one fasta file:

    cat ../one-to-one_orthogroups_Tfas_Tlei.perOG.txt | while read line ;
    do
     Orthogroup=`echo "$line"|awk '{print $1}'`;
     echo $Name;
     Tfas=`echo "$line"|awk '{print $2}'`;
     Tlei=`echo "$line"|awk '{print $3}'`;
     cat ../Tfas_seq/$Tfas:cds.fst ../Tlei_seq/$Tlei:cds.fst > $Orthogroup.fasta
    done

# Checking for completeness and length differences between orthologous pairs

Since Dn/Ds ratios are very sensitive to alignment quality, we first wanted to have an idea of the completeness of our gene pairs (i.e. presence of a start and stop codon) and their lengths. I therefore made a "checklist" containing this information for each gene using the script `script_make_checklist_completeness_length_orthologous_genes.py`

# Calculating Dn / Ds ratios for the final dataset
