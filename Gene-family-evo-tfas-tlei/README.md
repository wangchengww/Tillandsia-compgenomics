# Inferring Gene Family Evolution between *T. fasciculata* and *T. leiboldiana*

Using the inferred gene models of both reference genomes, we can get a first idea of which gene families have changed in size between the two species. However, before looking at the families with greatest changes, some filtering of faulty gene models has to be done first, especially for the *T. fasciculata* annotation.

# Detecting faulty gene models resulting from haplotigs in *T. fasciculata*

*T. fasciculata* is a wide-spread species, often regarded as a species complex, with very high levels of heterozygosity. Unsurprisingly, the accession used for our *T. fasciculata* reference genome was also more heterozygous than desired for genome assembly. This is especially an issue for inferring gene family sizes, since haplotig genes will be considered separate copies.
Using 50x whole genome sequencing data of our reference genome accession, we can look at coverage across the genome, and more importantly across gene models. The idea here is that legitimate duplications should have similar coverage to the rest of the genome. However, erroneous duplications that occurred for example due to the existence of haplotigs, will have a lower coverage than the rest of the genome, as reads from one gene will be distributed across two genic regions in the genome.

**1) Alignment of 50x WGS data to T.fas reference:**

Alignment of the raw data to T.fas genome was done using bowtie2. The choice of aligner is because the manual of bowtie2 explicitly mentions that reads aligning equally well in two different locations will be assigned to one of these two locations at random. This is important, since we want to infer a decrease in coverage when gene duplicates are false, and this would be more difficult if reads are mapped at multiple locations. BWA doesn't explicitly specify this.

    # Genome indexing:
    bowtie2-build -f Tillandsia_fasciculata_25_scaffolds.fasta Tillandsia_fasciculata_25_scaffolds
    # alignment
    bowtie2 --very-sensitive-local \
	-x /proj/grootcrego/Genome_assemblies/fasciculata/final/Tillandsia_fasciculata_25_scaffolds.fasta \
	-1 /proj/grootcrego/Genome_assemblies/fasciculata/0_raw_data/illumina/Tfas_Illumina_50x_trimmed_pair1.fq \
	-2 /proj/grootcrego/Genome_assemblies/fasciculata/0_raw_data/illumina/Tfas_Illumina_50x_trimmed_pair2.fq \
	-S Tfas_50x_illumina_to_Tfas25chrom.sam -p 24
    # Transforming to bam and sorting
    samtools view -Sb Tfas_50x_illumina_to_Tfas25chrom.sam > Tfas_50x_illumina_to_Tfas25chrom.bam
    samtools sort Tfas_50x_illumina_to_Tfas25chrom.bam -o Tfas_50x_illumina_to_Tfas25chrom.sorted.bam

**2) Calculating per-gene median and average coverage**

Per-base coverage in genic regions of T. fasciculata was calculated using samtools depth. I first obtained a bedfile with the coordinates of only our curated gene models (see orthofinder analysis):

    # Select only Tfas genes from file containing curated gene models
	grep "Tfasc_v1"  orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt > curated_Tfas_orthologues.txt
    # Obtain only gene ID
	cut -f 1 curated_Tfas_orthologues.txt > curated_Tfas_orthologues.IDonly.txt
    # Select curated genes from annotation GFF file
	grep -f curated_Tfas_orthologues.IDonly.txt /proj/grootcrego/Genome_assemblies/fasciculata/4_final_assembly/Tillandsia_fasciculata_v1.2.edited_allfeatures.gff > Tillandsia_fasciculata_v1.2.curated_orthologues_only.gff
    # Limit curated gff file only to mRNA entries (no exon, UTR, or CDS...)
	awk '$3 == "mRNA" {print $0}' Tillandsia_fasciculata_v1.2.curated_orthologues_only.gff > Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff
    # This curated GFF file was then converted to a BED format using gff2bed from bedops v. 2.4.37:  
    gff2bed < Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff > Tfas_ortholog_regions.gff.bed

Using this bedfile, I then calculated per-base coverage with samtools. Coverage of all bases included in the bedfile was reported, including areas of zero coverage (-a):

    samtools depth -a -b $bedfile $bamfile -o $output

I then appended the name of the gene at each position in the samtools depth output file using the python script `script_add_gene_names_to_cov_file.py`.

Note: This script was updated later on to include overlapping genes. See [here](https://github.com/cgrootcrego/Tillandsia-compgenomics/tree/main/Annotation/Gene_curation/overlapping_genes) for more info on these overlapping genes.

Mean and median coverage were then computed per gene and compiled into a table containing also orthology information with the script `script_calculate_mediancov_add_info.py`:

    python script_calculate_mediancov_add_ortho-info.py Tfas_genemodel_assessment/CovDepth_Tfas_ortholog_regions.fromgff.edited.txt Tfas_genemodel_assessment/curated_Tfas_orthologues.txt > computed_lengths_cov.txt

The file `computed_lengths_cov.txt` records the length of the vector of coverage entries for each gene. This was a way for me to make sure the script was recording the full length of overlapping genes.

# Correcting gene family sizes based on average coverage

I then made density plots of the average coverage per gene for different categories of genes using the Rscript `Assessing_multicopy_genemodels_cov.R`. These density plots showed that multi-copy gene models in *T.fasciculata* have a bimodal distribution of mean coverage - meaning that about half of the multicopy genes have a mean coverage around the genome-wide average, while the other half has a much lower mean coverage. I then determined a coverage threshold using a finite mixture model (FMM) with the package [cutoff].(http://marcchoisy.free.fr/fmm/index.html). This cutoff represents an estimated "split" of the two peaks in the bimodal distribution, in other words, the point at which, if we would split the bimodal distribution into two unimodal ones, a datapoint is unlikely to belong to both distributions. So, I used this threshold to separate "true" gene models from "faulty" ones. The threshold was determined at a median and mean coverage of 35. I then isolated all multicopy genes with a mean coverage under the threshold.

Note: the average full genome coverage was calculated by running samtools depth without specifying regions and running the awk one-liner `awk '{ total += $3; count++ } END { print total/count }`. To obtain the median coverage over the full genome, I ran `cut -f 3 CovDepth_perbase_fullgenome.txt | sort -n | awk -f median.awk`. Median.awk was the following short script:

    #/usr/bin/env awk
    {
    count[NR] = $1;
    }
    END {
    if (NR % 2) {
	    print count[(NR + 1) / 2];
    } else {
	    print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;
    }
    }

I then decided to rescale gene family sizes for the orthogroups which contained faulty multi-copy genes. There were 1981 such orthogroups, out of 17,641 (11 %), containing 6861 genes (26 %) I corrected gene family sizes by obtaining the total mean coverage in the orthogroup (sum of mean coverage of all Tfas genes in the OG) and dividing this by the "expected mean coverage" (full genome mean coverage x number of Tfas genes in the group). I called this ratio the correction factor, which I multiplied by the original number of   Tfas genes in the orthogroup to obtain a corrected family size.

For 63 orthogroups, the correction factor was > 1, meaning that the total average coverage of the orthogroup was larger than expected. In these cases, coverage is non-informative regarding the validity of genome sizes, so I decided not to correct for these families.

for 803 orthogroups, the total mean coverage of the orthogroup didn't even reach 46 (what we would expect for 1 gene). In this case, the correction usually brought the gene family size down to 1, but sometimes also to 0. In those cases, since the gene has been properly annotated and has an orthologous sequence in at least *T. leiboldiana* or *A.comosus*, I assumed the gene must exist and therefore corrected the size back up to 1.

I ran the same corrections for T. leiboldiana, although I couldn't separate faulty genes due to the fact the distribution was unimodal with just a small shoulder. However, while in the case of T.fas 1719 orthogroups were corrected, only 424 orthogroups were size corrected in T. leiboldiana.

Next, I integrated the new family sizes into my general per-gene orthology table, which we will feed back into gene family evolution analysis, with the python script `script_insert_corrected_sizes.py`.

# Studying gene family size differences in Tfas and Tlei

I then obtained a table of per-orhogroup gene counts by selecting the orthogroup and count fields of the orthology table.
The relationship of gene family size between species was explored in `Gene_family_evolution_new.R`

A first look at the most extreme families (high copy number in either species) showed that it may be interesting to exclude orthogroups containing chloroplastic / mitochondrial / ribosomal genes. To do that, I took the following steps:

    # Search for the words "chloroplast", "ribosomal" or "mitochondrial" and isolate orthogroups that carry such annotations
    grep -f filter_plastid_ribosomal_OGs.txt orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.txt | cut -f 7 | sort -u > mito_plastid_ribo_OGs.txt
    # Select all genes that don't belong to these orthogroups
    grep -v -f mito_plastid_ribo_OGs.txt orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.txt > orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.txt

Additionally, it also showed that in the case of leiboldiana, some corrections also have to be made. For example, OG15 counts 100 copies in Tlei, and seems a strong outlier. By checking his coverage, I came to the conclusion these copy numbers are probably also false.

# Changes to workflow July 27th

I made a few changes to the workflow above:
 - Size corrections were made for all multi-copy orthogroups in Tfas, instead of only on the orthogroups with genes < 34.5x. This only expanded the number of orthogroups run through size corrections from 1981 to 2632 (25 % more) and from 6861 genes to 8400 genes (19 % more).
 - Size corrections were also made for upward corrections in orthogroups that did not belong to chloroplasts.
 - Size corrections were made not based on the full genome average but on the average coverage of ancestral single-copy genes. The idea is that this category represents better the sort of coverage expected in genic regions (full genome coverage is extremely variable). Additionally, I used the 25th and 75th percentile of ancestral single-copy genes to devise an interval within which corrections were not made to account for variability in coverage.
