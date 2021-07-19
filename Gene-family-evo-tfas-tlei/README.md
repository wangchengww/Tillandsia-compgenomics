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

I then made density plots of the average coverage per gene for different categories of genes using the Rscript `Assessing_multicopy_genemodels_cov.R`. These density plots showed that multi-copy gene models in *T.fasciculata* have a bimodal distribution of median coverage - meaning that about half of the multicopy genes have a median coverage around the genome-wide average, while the other half has a much lower median coverage. I then determined a coverage threshold using a finite mixture model (FMM) with the package [cutoff].(http://marcchoisy.free.fr/fmm/index.html). This cutoff represents an estimated "split" of the two peaks in the bimodal distribution, in other words, the point at which, if we would split the bimodal distribution into two unimodal ones, a datapoint is unlikely to belong to both distributions. So, I used this threshold to separate "true" gene models from "faulty" ones. The threshold was determined at a median coverage of 35. I then isolated all multicopy genes with a median coverage under the threshold.

Note: the average full genome coverage was calculated by running samtools depth without specifying regions and running the awk one-liner `awk '{ total += $3; count++ } END { print total/count }`. To obtain the median coverage over the full genome, I ran `sort -n CovDepth_perbase_fullgenome.txt | awk -f median.awk`. Median.awk was the following short script:

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

I then decided to rescale gene familz sizes for the orthogroups which contained faulty multi-copy genes. There were 6927 such orthogroups. I corrected gene family sizes by obtaining the total median coverage in the orthogroup (sum of median coverage of all Tfas genes in the OG) and dividing this by the "expected median coverage" (full genome median coverage x number of Tfas genes in the group). 
