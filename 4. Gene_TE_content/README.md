# Visualizing genic and repetitive content in *T. fasciculata* and *T. leiboldiana* genomes

Here are all the steps to produce the circular figure showing window-based repetitive and genic content for each scaffold of both species' genome assemblies.

# Calculating genic content

For this, the final orthology file orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt (see Orthology) was used. Gene content per 1 MB window was calculated both as the number of genes starting per window and the proportion of bases being in genic regions. This was done for all orthologous genes, only 1-1 genes and only multicopy genes for both species in the Rscript `Calculate_Gene_content.R`.

# Calculating repetitive content

Using the masked version of the assemblies (run after extensive EDTA run, see Synteny, LastZ runs), I calculated the percentage of masked bases per 1 MB windows with the python script `script_calculate_repetitive_content_perwindow_maskedfasta.py`:

    python calculate_repetitive_content_perwindow_maskedfasta.py Tillandsia_leiboldiana_26_scaffolds.fasta.masked 1000000

# Visualizing per-scaffold genic and repetitive content

The circular figures showing per-scaffold genic and repetitive content were made using circlize in R with the script `Visualize_Genic_TE_content_circlize.R`.

# In depth analysis of repetitive content

See `Repetitve_content_in_depth_analysis.md`.

# Calculating per scaffold genic / repetitive content ratios

Since visually it is clear that TE content is larger in T. leiboldiana, I decided to compute genic-to-repetitive content ratios for each contig for both species, to have a numeric comparison. We already have the number of repetitive basepairs per scaffold from the in-depth analyses of TE content (see above). To obtain total numbers of genic basepairs per scaffold, I counted basepairs marked as exon from the gff files of both genomes with the script `script_genic_proportion_perscaff.py`. Importantly, these were run on subsets of the gff file containing only curated orthologues: `python ../Tfas_assembly/script_genic_proportion_per_scaff.py Tillandsia_leiboldiana_v1.2.curated_orthologues_only.gff chrNameLength.txt`
