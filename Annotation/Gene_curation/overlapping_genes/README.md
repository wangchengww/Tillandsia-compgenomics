# Investigating and curating overlapping gene models

When studying gene families in *T. fasciculata*, I realized that a number of annotated gene models had overlapping ranges in the genome. A bit of research unveiled that this commonly happens in MAKER when several gene prediction algorithms are applied (see [here](https://groups.google.com/g/maker-devel/c/I3qYzziurso)). This often leads to overlapping gene models being called on opposite strands.

 I decided to check more in depth how many genes on the plus and minus strand overlap with each other and how we could deal with this issue. I only worked with our curated set of genes for this:

    awk '$7 == "+" {print $0}' Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff > Tfas_curated_orthologs_mRNA_PLUS.gff
	awk '$7 == "-" {print $0}' Tillandsia_fasciculata_v1.2.curated_orthologues_only.mRNA.gff > Tfas_curated_orthologs_mRNA_MINUS.gff
	bedtools intersect -a Tfas_curated_orthologs_mRNA_PLUS.gff -b Tfas_curated_orthologs_mRNA_MINUS.gff -wo > Tfas_overlapping_genemodels_opposite_strands.txt

This revealed 705 pairs of overlapping genes. The bedtools intersect output file is not very easy to navigate, so I reformatted it with the python script `reformat_overlap_file.py`.
When navigating this table, it becomes clear that the 705 overlapping gene pairs belong to all kinds of categories (single-copy, multi-copy, unique to T.fas ...) and that there is most of the time some form of evidence for each gene (see QI and eAED scores). The 705 genes on the PLUS strand, belong to 601 different orthogroups, on the MINUS strand 605 orthogroups are represented. I calculated if overlapping gene pairs were from the same orthogroup or not with a short bash loop:

    i=0
	j=0
	cat Tfas_overlapping_genemodels_opposite_strands.table.txt | ( while read line; do
	og1="$(echo "$line"|awk 'BEGIN { FS="\t" } {print $2}')"
	og2="$(echo "$line"|awk 'BEGIN { FS="\t" } {print $12}')"
	if [[ "$og1" == "$og2" ]] ;
	then
		i=$((i+1))
	else
		j=$((j+1))
	fi
	done
	echo "$i gene pairs are from the same orthogroup"
	echo "$j gene pairs are from different orthogroups" )

This revealed that only 4 gene pairs were from the same OG, all remaining 701 are from different OGs. Additionally I checked how many gene pairs had a gene fully overlapping with the other gene with the following one liner:
    awk 'BEGIN { FS="\t" } $6 < $20 || $15 < $20 {print $0}' Tfas_overlapping_genemodels_opposite_strands.table.txt | wc -l

This resulted in 245 fully overlapping genes (35 %).
Overall, it makes it difficult to curate this group of genes, since there is no clear signal that any of these genes should be removed. Overlaps are often only a few hundred basepairs, so it could also be an imprecise gene boundary by MAKER. For fully overlapping genes, the AED score is still relatively high. So we either have to fully remove these genes or keep them.
