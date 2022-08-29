# Inferring signatures of selection by calculating dN/dS ratios between T. fasciulata and T. leiboldiana

In this section, we infer dN / dS ratios for all single-copy pairs of orthologous genes between T. fasciculata and T. leiboldiana. Part of this work was done by Francesca Beclin. We obtained Dn / Ds ratios by aligning the orthologous sequences using a codon-aware aligner (MACSE) and inferring ratios with KaKs_Calculator and PAML.

# Selecting single-copy orthogroups and compiling per-orthogroup fastafiles

First I compiled the sequences of each gene into per-orthogroup fastafiles. For this, I extracted the one-to-one genes obtained in orthofinder analyses from the annotation gff files of each assembly. The orthofinder file used for this was: /scratch/grootcrego/orthofinder/run_orthofinder_Tfas_Tlei_Acom_25_scaffolds/OrthoFinder/Results_Jan22/Phylogenetic_Hierarchical_Orthogroups/one-to-one_orthogroups_Tfas_Tlei.txt and was created in the postprocessing of orthofinder results (see Orthology&Synteny notes). There are 13,128 orthologous pairs in this file. This file is a subset of orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt

To extract the features belonging only to one-one orthologous genes:

    grep -f one-to-one_orthogroups_Tfas_Tlei.IDonly.txt \
    Tillandsia_fasciculata_v1.2.edited_allfeatures.gff > \
    Tillandsia_fasciculata_v1.2.edited_allfeatures.one-to-one_orthologs.gff

Then I extracted the coding sequence features and replaced ID with Name to make the below python script work:

    awk '$3 == "CDS" {print $0}' [gff] | sed 's/ID=/Name=/g' > [gff.CDSonly]

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
This file was eventually modified in R to a per-orthogroup format in which the absolute difference in sequence length between pairs was calculated. Additionally, Francesca calculated relative differences for each species (diff/length_Tfas and diff/length_Tlei) and noted completeness in the orthogroup (both complete / one complete / zero complete). Additionally, francesca counted the exon number per gene and the difference in counts between genes.

# Testing filtering options of alignments based on length differences and completeness

We created four subsets with different degrees of filtering to run preliminary dN/dS analyses on:

- stringent filtering: both alignments are complete (they have a start and stop codon), they have an absolute length difference of maximum 200 basepairs and they have a relative length difference that is smaller than 0.1 (meaning that the difference in length is no more than 10 % of the total length of that gene).
- relaxed completeness: same requirements as above regarding length differences, but only one gene is complete.
- relaxed relative difference: both genes have to be complete and the same requirement applies for absolute length differences as above. However, the relative length difference ranges from 0.1 to 0.2.
- relaxed completeness and relative difference: absolute length requirement remains the same, but relative differences range from 0.1 to 0.2 and only one gene is complete.

The total number of orthogroups for each subset was:
6208  stringent filtering (47 %)
493   relaxed relative difference (4 %)
2569  relaxed completeness (19 %)
316   relaxed relative difference and completeness (2 %)

Francesca then randomly selected 200 orthogroups for each filtering setting and performed alignments with them.

Alignments on subsets were performed twice, once with very basic code and once with more specific requirements that optimize the alignment:

    # Basic alignment
    for i in *; do
     java -jar macse_v2.05.jar -prog alignSequences -seq $i ;
    done
    # Extensive alignment
    for i in *; do
     java -jar macse_v2.05.jar -prog alignSequences -seq $i -local_realign_init 1 -local_realign_dec 1 ;
    done

The alignments were then converted into .axt files using a modified form of the script
03_ConvertFasta.py from the {AlignmentProcessor}(https://github.com/WilsonSayresLab/AlignmentProcessor/blob/master/bin/03_ConvertFasta.py) software. Preliminary dN/dS ratios were then calculated with KaKs calculator and the MYN model using a modified script of 04_CallKaKs.py from ALignmentProcessor.

Comparing the distribution of dN/dS ratios between filtering categories showed that relaxing length difference didn't cause extreme values of dN/dS. However, relaxing completeness did cause a few estimates to show up, though most estimates were still in the same range as with stringent filtering. I therefore decided to keep the filtering loose at this point, and only remove alignments of genes where both are incomplete and where length differences are > 0.2 in both genes. This kept 10,362 orthologous pairs for pairwise alignment, which is almost 80 % of the original set.

# Calculating pairwise dN / dS ratios for the final dataset

Next, using the checklist, I selected all genes:

	grep -v "0 complete" checklist_with_diff_completeness.txt | \
	awk '$5 > -0.2 && $5 < 0.2 && $6 > -0.2 && $6 < 0.2' | cut -f 1 | \
	sed 's/"//g' > orthogroups_final_subset_0.2_lengthdiff_1_complete.txt

Alignment was done with optimization parameters for all orthologous pairs with the bash script `run_macse.sh`.

As CodeML allows for pairwise calculation of dN / dS values, I decided to drop KaKsCalculator and find a way to automate CodeML. The selected alignments were converted to PHYLIP format using `03_ConvertFastatoAXTorPhylip.py`. I used modified scripts from ALignmentProcessor `04_CallCodeML_modified.py` and `parallelcodeML.py` to run codeML automatically for all 10,362 orthologous pairs.

codeML was run twice for a null model and alternative hypothesis. The settings are summarized in the two control files `codeml_null.ctl` and `codeml_alt.ctl`. Concretely, we used the site-specific model MO (Nssties = 0) for both runs. In the null model, omega (dN/dS) was fixed to 1, whereas in the alternative, it was estimated from the data starting from 1. The reason to run codeML twice is that we gather likelihoods for the null and alternative model and can compute p-values, to see how significant our inferences are.
Additional choices in the settings of codeML were to set codon frequencies to be inferred from the data at all three positions in the codon (F3x4) (CodonFreq = 2) and kappa (Ts/Tv ratio) was set to 3, but is estimated from data (kappa = 3, fix_kappa = 0). The results of both runs were compiled with the script `script_compile_codeml_LRT.py`. This script also performs a likelihood ratio test (chisquare). I called the script on all files as follows:

    null=(/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/dnds/dnds_calculations/codeML/null/*.out)
    alt=(/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/dnds/dnds_calculations/codeML/alt/*.out)
    n=0
    for ((i=0;i<=${#null[@]};i++))
    do
      python script_compile_codeml_LRT.py "${null[i]}" "${alt[i]}"
      n=$((n+1))
      echo $n
    done

Important: I modified the phylip files slightly, because codeML doesn't accept "!" in the alignments. MACSE introduces these whenever there is a change in frameshift. In other words, when an entire codon is deleted, this will be shown as "---" but when there is a gap < 3, it will show as "!!A" or something of the like. I replaced all "!" by "-" so that codeML wouldn't throw errors.

The dN/dS results were then analyzed with `Assessment_dnds_values.R`. Alignments of significant genes were checked with AliView. Additionally, RNA-seq data used for annotation was mapped back to the main scaffolds and visualized in Jbrowse to further assess alignments.

# Curating candidate genes

One danger of dN/dS is that this measure is highly dependent on alignment quality and high values of dN/dS can often be due to misalignment. I had a manual, visual look at all ~ 40 candidate alignments and detected many odd alignments. Since solving alignments can be tricky on a visual basis, I decided to rerun alignments with a different pipeline, by aligning the protein sequences with ClustalOmega and then running PAL2NAL (which will automatically compute dN/dS in PAML).

First, I extracted the gene names per orthogroup for all 38 groups from the orthogroup file. Then, I created protein fasta files per orthogroup:

	cat list_Candidate_genes_per_OG.txt | while read line; do
		Orthogroup=`echo "$line"|awk '{print $3}'`;
		echo $Orthogroup;
		Tfas=`echo "$line"|awk '{print $1}'`;
		Tlei=`echo "$line"|awk '{print $2}'`;
		seqkit grep -p $Tfas /gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/Tfas_assembly/Tillandsia_fasciculata_v1.2.protein.fasta > $Orthogroup.AA.fasta
		seqkit grep -p $Tlei /gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/Tlei_assembly/Tillandsia_leiboldiana_v1.2.protein.fasta >> $Orthogroup.AA.fasta
   done

Next, I made a pairwise protein alignment with clustalo:

	for i in ../protein_seq/*; do
		name=$(basename "$i" .AA.fasta)
		echo $name
		clustalo -i $i -o ./$name.alignment -t Protein --outfmt clu
	done

This was then used in pal2nal to obtain a final nucleotide alignment:

	cat ../list_Candidate_genes_per_OG.txt | while read line; do
		Orthogroup=`echo "$line"|awk '{print $3}'`;
		echo $Orthogroup;
		pal2nal.pl ../protein_alignments/$Orthogroup.alignment ../nuc_fasta/$Orthogroup.fasta -nogap -output paml > $Orthogroup.pal2nal.out
	done

codeML was then run for the null and alternative models as explained above, and a summary file with Likelihood ratio test was compiled as well with python script. Out of the 38 candidate orthogroups detected using MACSE alignments, only 19 were kept using pal2nal. After seeing these results, I decided to run pal2nal for all genes.

# Rerunning alignments with PAL2NAL

This was done as specified above, and resulted in 44 gene pairs with significant dN/dS > 1. 19 (44 %) of these overlap with MACSE. Since I ran pal2nal on the entire set of 1-1 orthologues (without any filtering), I also obtained 16 candidate genes that should've been filtered out. I guess this reinforces the need to filter to avoid false positives. Additionally, 9 alignments came out as candidates from pal2nal alignment that were not significant with MACSE. Here, it mostly seems due to uncertainties in how to align gaps. I therefore decided to only work with the 19 pairs that were called with both alignment methods.

# Obtaining Tillandsia outgroup sequences for branch-site models

Next, we decided to recalculate dN/dS values using branch-site models. This allows omega to vary across branches and can give us an idea of which species experienced adaptive sequence evolution in each gene, as pairwise calculations cannot give us direction. However, this means the inclusion for an outgroup. While Ananas comosus gene models are readily available from our orthofinder analyses, this is a widely diverged taxon compared to our Tillandsias and may cause alignment issues. Therefore, we decided to make use of existing RNA-seq of an outgroup Tillandsia to obtain outgroup sequences.

Data availability:
- Tillandsia floribunda: CAM plant of the subgenus Allardtia, part of the T. biflora complex. 2 RNA-seq samples are available.
- Tillandsia australis: C3 plant of the subgenus Allardtia. 3 RNA-seq samples are available.
- Tillandsia sphaerocephala: CAM plant of the subgenus Allardtia, part of the T. biflora complex. 3 RNA-seq samples are available.

Since all three species are roughly equally distant to our Tillandsias, there is no obvious favourite to include as an outgroup. It might be interesting to include two outgroups, where one is a C3 and another a CAM species (T. australis and T. sphaerocephala). I eventually decided to work with T. australis and T. sphaerocephala, as we have more data available.

The RNA-seq data of day and night sample of each individual were combined into one file: `cat Taus_C3_day_A_R2.fastq Taus_C3_night_A_R2.fastq > Taus_C3_A.combined.R2.fastq`

Trimming was run with AdapterRemoval with following loop:

    R1=($input_dir/*.R1.fastq)
    R2=($input_dir/*.R2.fastq)
    for ((i=0;i<=${#R1[@]};i++))
    do
      AdapterRemoval --file1 "${R1[i]}" --file2 "${R2[i]}" --basename "${R1[i]%.R1.fastq}" --trimns --trimqualities --minquality 20 --trimwindows 12 --minlength 36
      fastqc "${R1[i]%.R1.fastq}.pair1.truncated" "${R2[i]%.R2.fastq}.pair2.truncated"
    done

After trimming, samples contained between 33.2 - 40.9 million pairs.

Mapping was then performed against the T. leiboldiana genome (fewer multimappers) with STAR with the following loop:

    # Index Genome
    $star --runMode genomeGenerate --genomeDir $ref_dir --genomeFastaFiles $ref --sjdbGTFfile $gff_file --sjdbGTFtagExonParentTranscript Parent  --sjdbOverhang 149 --runThreadN 24
    # Map samples to reference
    R1=($wd/*.pair1.truncated)
    R2=($wd/*.pair2.truncated)
    for ((i=0;i<=${#R1[@]};i++))
    do
     /apps/star/2.5.3a/bin/Linux_x86_64/STAR --genomeDir $ref_dir \
     --readFilesIn "${R1[i]}" "${R2[i]}" --outSAMtype BAM Unsorted \
     --outFileNamePrefix "${R1[i].%.pair1.truncated}" \
     --limitBAMsortRAM 10240000000 --runThreadN 48
     done

Unique mapping rates ranged between 56.9 % to 72.9 %.

Tibo then took over the work to extract fasta sequences for candidate genes from all individuals of T.aus and T.sphae. This resulted in a per-gene multi fasta containing two fastasequences for each individual. The second sequence contains non-ref heterozygous (unfixed) variants. Therefore, we will always work with the first sequence as dN/dS rates are based on fixed variants.

To choose which individual to use for each gene, I counted the number of N's: `seqtk comp  [gene.fasta] | cut -f 1,9`. This is because Tibo noted some complications leading to increased N's in the sequences and therefore the fewer N's the better, especially because length tends to be the same.

For the test, I decided to work with orthogroup OG0006131. Here, in T. australis, individual C has the fewest N's (7), and for T.sphaerocephala, the number of N's is equal in all individuals so I decided to work with A as it has the highest mapping rates.

This part of the work was eventually dropped.

# Running dN/dS for 1:1:2 and 1:2:1 paralogs

I decided to also run parallel pairwise tests for duplicated orthogroups that are either 1:1:2 (duplicated in Tlei) or 1:2:1 (duplicated in Tfas).

There are 193 orthogroups in 1:1:2 conformation, and 917 orthogroups in 1:2:1 conformation. However, when looking at the corrected orthogroup sizes, we have 175 groups in 1:1:2 and 263 in 1:2:1. I decided to be very stringent with this analysis and only work with orthogroups of the intersection of corrected and uncorrected sizes (the group is reported as 1:1:2 or 1:2:1 both in the corrected and uncorrected files). This resulted in 108 groups in 1:1:2 and 190 groups in 1:2:1.

Obtaining fasta-sequences per orthogroup was done similarly as above but with a few changes. For each orthogroup, we end up with two fastafiles, for example in orthogroups of 1:1:2, we get one file with the Tfas gene and copy 1 of the Tlei gene, and another file with the Tfas gene and copy 2 of the Tlei gene. For this, I had to modify the python script to include three columns, into `script_compile_per-orthogroup_gene_list.paralog.py`. Concatenation of the fastasequences into two files happened with following loop:

	cat ../orthologs-121.perOG.txt | while read line ;
	do
 		Orthogroup=`echo "$line"|awk '{print $1}'`;
 		echo $Orthogroup;
 		Tfas1=`echo "$line"|awk '{print $2}'`;
 		Tfas2=`echo "$line"|awk '{print $3}'`;
 		Tlei=`echo "$line"|awk '{print $4}'`;
 		cat ../Tfas_seq/$Tfas1:cds.fst ../Tlei_seq/$Tlei:cds.fst > ${Orthogroup}_Tfas-copy1.fasta
 		cat ../Tfas_seq/$Tfas2:cds.fst ../Tlei_seq/$Tlei:cds.fst > ${Orthogroup}_Tfas-copy2.fasta
	done

Because we start with a small number of orthogroups, I decided not to do any previous filtering. So I ran macse with a modified version of run_macse.sh for 1:1:2 and 1:2:1 genes.

Alignments were converted to phylip format with the above modified python script:

    python 03_ConvertFastatoAXTorphylip_modfied.py -i alignments_121/ -o phylip_121/ --phylip  

And run as above with a null and alt model in codeML. I compiled the results and performed LRT with a modified version of the above python script, `script_compile_codeml_LRT.paralogs112.py`. This revealed only one orthogroup with a significant dn/ds >1 in 1:1:2 and two orthogroups in 1:2:1. These are summarized in the google spread sheet. Unfortunately, the XAP circadian timekeeper which exists in 1:2:1 conformation did not show a dN/dS > 1, which would've been interesting in combination with the RNA-seq results, where we see that different copies are recruited in either species.

# dN/dS of Phosphoenolpyruvate carboxylase (PEPC)

PEPC occurs in 2:2:2 conformation and only shows significant expression in both species in one of the two copies. Therefore, one could hypothesize that the other copy is not functional anymore. I thought that I could run dN/dS calculations between all four gene sequences (2 for Tlei and 2 for Tfas) to see if we see major deviation from dN/dS, especially between the highly expressed copy in Tfas and the two copies in Tlei.

# Obtaining per-chromosome dN/dS statistics

Apart from candidate genes, I also wanted to see the general dN/dS distribution across chromosomes. To look at this, I decided to only work with alignments that had at least 5 variant sites. I obtained this information with AMAS and created a list of orthogroups to keep. This resulted in 9,077 orthogroups.
