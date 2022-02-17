# dN / dS Branch-specific estimates for candidate genes

# Files in folder

1. *.bam : RNA-seq data of three accessions of T. australis and T. sphaerocephala, aligned to the T. leiboldiana genome. These are C3 ad CAM species respectively of a Tillandsia clade related to the subgenus.

2. Tillandsia_leiboldiana_26_scaffolds.fasta: Unmasked reference genome containing only the 26 main scaffolds

3. Tillandsia_leiboldiana_v1.2.edited_allfeatures.26chrom.gff: Full GFF file of the T. leiboldiana reference genome including functional information

4. Genes_one-to-one_dNdS_subset.txt: List of 20,724 genes previously used in pairwise dN/dS calculations.

5. Candidate_genes_dNdS.txt: 20 genes that had dN/dS > 1 and p < 0.05 in pairwise tests.

# To do:

1. Call VCF from BAM Files
2. Extract genome sequence
3. Cut out sequences of candidate genes
4. Align outgroup sequences and gene models
5. Run branch-specific model to estimate dN/dS
