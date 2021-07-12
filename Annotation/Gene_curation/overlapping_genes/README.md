# Investigating and curating overlapping gene models

When studying gene families in *T. fasciculata*, I realized that a number of annotated gene models had overlapping ranges in the genome. A bit of research unveiled that this commonly happens in MAKER when several gene prediction algorithms are applied (see [here](https://groups.google.com/g/maker-devel/c/I3qYzziurso)). This often leads to overlapping gene models being called on opposite strands.

 I decided to check more in depth how many genes on the plus and minus strand overlap with each other and how we could deal with this issue. I only worked with our curated set of genes for this:

    awk '$7 == "+" {print $0}' Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff > Tfas_curated_orthologs_mRNA_PLUS.gff
	awk '$7 == "-" {print $0}' Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff > Tfas_curated_orthologs_mRNA_MINUS.gff
	bedtools intersect -a Tfas_curated_orthologs_mRNA_PLUS.gff -b Tfas_curated_orthologs_mRNA_MINUS.gff -wo > Tfas_overlapping_genemodels_opposite_strands.txt

This revealed 705 overlapping genes. The bedtools intersect output file is not very easy to navigate, so I reformatted it with the python script `reformat_overlap_file.py`.
