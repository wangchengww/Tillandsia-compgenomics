# Calling structural variants and comparing caller softwares

To study gene family evolution in the subgenus Tillandsia, one possible strategy would be to call family sizes in several Tillandsia species, including *T. leiboldiana* and *T. fasciculata*. In order to do this, we would have to call structural variants, which include deletions and duplications, and find variants that overlap with genes.

However, our available data mostly consists of 10x-20x illumina data, and I first wanted to see how reliable structural variant callers were using such data. I therefore made a comparative analysis between two PacBio SV callers and three illumina-based callers using the same reference (T. leiboldiana) and the same sample (T. fasciculata), for which we luckily had both types of data.

Our Pacbio data was about 35x coverage, and illumina data was 50x coverage. I also wanted to see the effect of different coverage levels, so structural variants were called with illumina data at 5x, 10x, 20x and 50x.

# Calling structural variants with SVIM

SVIM is able to call interspersed duplications, making it quite an interesting software. I ran SVIM with the following command:

`svim reads Tfas_Pacbio_to_Tlei DTG-DNA-541.subreads.fasta Tillandsia_leiboldiana_26_scaffolds.fasta`

# Calling structural variants with Sniffles

Sniffles was run to have some form of control within the PacBio methods. Though we assume PacBio SV calling will be more accurate, I realized after an initial run that some of the output of SVIM was a bit irregular (overlapping deletions), so I included a second software:

`pacbio_bam=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/calling_SV/sniffles/DTG-DNA-541.subreads.ngmlr.coordsorted.bam
output=StrucVar_Tfas_to_Tlei_sniffles.vcf
wd=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/calling_SV/SVIM

cd $wd
sniffles -t 48 -m $pacbio_bam -v $output`

# Subsetting coverage of illumina fastq files

Using rasusa, I randomly selected reads to ensure an average 5x, 10x, 20x coverage for each file. The original data was sequenced at 50x.

First I estimated the genome size:  
`grep -v ">" Tillandsia_leiboldiana_26_scaffolds.fasta | wc | awk '{print $3-$1}'`

This resulted in 906,467,929 bp. Subsets based on coverage were then created like so:  
`rasusa --coverage 20x --genome-size 906467929 --input Tfas_Illumina_50x_trimmed_pair1.fq Tfas_Illumina_50x_trimmed_pair2.fq --output Tfas_Illumina_50x_trimmed_pair1.fq Tfas_Illumina_50x_trimmed_pair2.fq`

# Calling SV with Delly

Preparation:
I prep the input files by aligning them, add read groups, marking duplicates, sorting and indexing.

`bwa mem -t 8 $ref_genome $pair1 $pair2 | samtools view -Sb - | samtools sort -@4 - -o Tfas_illumina_to_Tlei_ref.10x.sorted.bam  
for file in $filedir ; do  
 	filename="$(basename $file)"  
	echo $filename  
	picard AddOrReplaceReadGroups I=$file o="$(basename $file)".RG.bam RGLB=WGD RGPL=illumina RGPU=Lib1 RGSM=T.fasciculata_B1840 RGID=T.fasciculata_B1840  
	picard MarkDuplicates  I=$filename o=${filename%.bam}.Dup.bam M=${filename%.bam}.Dup.bam_metrics.txt  
	samtools sort -o ${filename%.bam}.Dup.sorted.bam ${filename%.bam}.Dup.bam  
	samtools index ${filename%.bam}.Dup.sorted.bam  
done`  

Running Delly:

`$delly call -g $ref_genome -o $output.bcf $bam`
