# Visualizing genic and repetitive content in *T. fasciculata* and *T. leiboldiana* genomes

Here are all the steps to produce the circular figure showing window-based repetitive and genic content for each scaffold of both species' genome assemblies.

# Calculating genic content

For this, the final orthology file orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt (see Orthology) was used. Gene content per 1 MB window was calculated both as the number of genes starting per window and the proportion of bases being in genic regions. This was done for all orthologous genes, only 1-1 genes and only multicopy genes for both species in the Rscript `Calculate_Gene_content.R`.

# Calculating repetitive content

Using the masked version of the assemblies (run after extensive EDTA run, see Synteny, LastZ runs), I calculated the percentage of masked bases per 1 MB windows with the python script `script_calculate_repetitive_content_perwindow_maskedfasta.py`:

    python calculate_repetitive_content_perwindow_maskedfasta.py Tillandsia_leiboldiana_26_scaffolds.fasta.masked 100000
