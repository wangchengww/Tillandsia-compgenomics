# Annotation of repetitive content in *T. fasciculata* and *T. leiboldiana*

TE annotation was performed using the Extensive *de novo* TE Annotator ({EDTA}(https://github.com/oushujun/EDTA)) for both species. This was after Gil performed a few tests comparing RepeatModeller and EDTA showing that BUSCO scores were higher when using EDTA.

# Original annotation, used for masking during gene prediction

Originally, EDTA was run in version 1.7.8 by Gil for Tillandsia fasciculata (on the cube) and version 1.9.5 by me for Tillandsia leiboldiana (on the VSC4). The original annotation didn't include extensive corrections as these were not available in those versions. Annotation was run with the command:

    perl EDTA.pl --genome tillandsia_leiboldiana_assembly.pilon.upper.sorted.fasta \
	--cds pineapple.20150427.cds.wo_transposons.fasta \
	--sensitive --threads 16

*A. comosus* coding sequences were used to filter out genes falsely identified as TEs in the final step of annotation.

Important: The original annotation files for T. fasciculata were lost during a clean-up of the cube. However, the library (still available) was used in the downstream gene annotation process as well as masked versions of the genome (still available).

# Detailed annotation of the full assembly

To report the repetitive content of both assemblies, I reran annotation with extra steps to obtain the most precise result. This was done on the cube for both full assemblies as the VSC4 was throwing complex errors. The EDTA version used here was 1.9.6:

	perl $edta --genome $genome --cds $cds --step all \
	--sensitive 1 --anno 1 --evaluate 1 --threads 24
