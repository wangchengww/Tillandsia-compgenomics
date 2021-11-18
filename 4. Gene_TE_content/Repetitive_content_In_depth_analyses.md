# In depth analyses of repetitive content in *T. fasciculata* and *T. leiboldiana*

# Calculating abundances of TE classes for the whole genome and per scaffold

Since I was facing issues with EDTA running to completion, I devised the python script `script_TE_abundances.py` which will count the number of basepairs per scaffold and for the whole genome of each TE class that EDTA identifies.

Because I eventually did manage to run the --anno step in EDTA, I had the .TEanno.sum file to compare my script to. The percentages for the whole genome are similar but can differ by about 1 percentage point due to how I deal with nested TEs. Nested TEs were not counted in the basepair count of each class, since this would make the script very complicated. This explains probably the slight difference between my method and EDTA's method

I assessed per-scaffold TE content in the Rscript `Analyze_per-scaffold_TE_abundances.R`. This showed that while large scaffolds had on average 65 % TE content in T. fasciculata and ... in T. leiboldiana, small scaffolds consisted on average of 94 % TEs. 35 % of all TE content can be found on small scaffolds, which make up 25 % of the assembly. This, together with evidence from orthofinder showing that almost all orthologous genes are on the main scaffolds, shows that small scaffolds in the assembly are mostly TE remnants that were not incorporated in the main assembly. 
