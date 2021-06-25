# Calling structural variants and comparing caller softwares

To study gene family evolution in the subgenus Tillandsia, one possible strategy would be to call family sizes in several Tillandsia species, including *T. leiboldiana* and *T. fasciculata*. In order to do this, we would have to call structural variants, which include deletions and duplications, and find variants that overlap with genes.

However, our available data mostly consists of 10x-20x illumina data, and I first wanted to see how reliable structural variant callers were using such data. I therefore made a comparative analysis between two PacBio SV callers and three illumina-based callers using the same reference (T. leiboldiana) and the same sample (T. fasciculata), for which we luckily had both types of data.

Our Pacbio data was about 35x coverage, and illumina data was 50x coverage. I also wanted to see the effect of different coverage levels, so structural variants were called with illumina data at 5x, 10x, 20x and 50x.

# Calling structural variants with SVIM

SVIM is able to call interspersed duplications, making it quite an interesting software. I ran SVIM with the following command:

`svim reads Tfas_Pacbio_to_Tlei DTG-DNA-541.subreads.fasta Tillandsia_leiboldiana_26_scaffolds.fasta`
