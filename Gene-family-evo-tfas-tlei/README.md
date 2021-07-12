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

We lose 3 genes here out of 26,325 because they somehow have overlapping regions with other genes. The missed genes were: Tfasc_v1.24295-RA, Tfasc_v1.02150-RA and Tfasc_v1.03676-RA.
I realized the genes overlap with other genes on opposite strands, and by doing some research saw this is relatively common in maker annotations when different predictors are implemented. Therefore, I decided to check how many genes on the plus and minus strand overlap with each other:

    awk '$7 == "+" {print $0}' Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff > Tfas_curated_orthologs_mRNA_PLUS.gff
	awk '$7 == "-" {print $0}' Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff > Tfas_curated_orthologs_mRNA_MINUS.gff
	bedtools intersect -a Tfas_curated_orthologs_mRNA_PLUS.gff -b Tfas_curated_orthologs_mRNA_MINUS.gff -wo

This revealed 705 overlapping genes. I have to look into these genes and decide what to do with them in the future.
