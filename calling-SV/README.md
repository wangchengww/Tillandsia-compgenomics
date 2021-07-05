# Calling structural variants and comparing caller softwares

To study gene family evolution in the subgenus *Tillandsia*, one possible strategy would be to call family sizes in several  species, including *T. leiboldiana* and *T. fasciculata*. In order to do this, we would have to call structural variants, which include deletions and duplications, and find variants that overlap with genes.

However, our available data mostly consists of 10x-20x illumina data, and I first wanted to see how reliable structural variant callers were using such data. I therefore made a comparative analysis between two PacBio SV callers and three illumina-based callers using the same reference (*T. leiboldiana*) and the same sample (*T. fasciculata*), for which we luckily had both types of data.

Our Pacbio data was about 35x coverage, and illumina data was 50x coverage. I also wanted to see the effect of different coverage levels in illumina data, so structural variants were called with illumina data at 5x, 10x, 20x and 50x.

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

`bwa mem -t 8 $ref_genome $pair1 $pair2 | samtools view -Sb - | samtools sort -@4 - -o Tfas_illumina_to_Tlei_ref.10x.sorted.bam \
for file in $filedir ; do \
 	filename="$(basename $file)" \
	echo $filename \
	picard AddOrReplaceReadGroups I=$file o="$(basename $file)".RG.bam RGLB=WGD RGPL=illumina RGPU=Lib1 RGSM=T.fasciculata_B1840 RGID=T.fasciculata_B1840 \
	picard MarkDuplicates  I=$filename o=${filename%.bam}.Dup.bam M=${filename%.bam}.Dup.bam_metrics.txt \
	samtools sort -o ${filename%.bam}.Dup.sorted.bam ${filename%.bam}.Dup.bam \
	samtools index ${filename%.bam}.Dup.sorted.bam \
done`  

Running Delly:

`$delly call -g $ref_genome -o $output.bcf $bam`

# Calling SV with lumpy

Preparation:
First, I aligned the separate fastq files, added read groups and removed duplicates as recommended by the manial. The discordant and split-read alignments were extracted. Then this was all sorted and indexed.

`for ((i=0;i<=${#pair1[@]};i++))
do
  Name=$(basename ${pair1[i]})
  Name=${Name%_trimmed_pair1.fq}
  echo $Name
  bwa mem -t 48 -R "@RG\tID:$sample_id\tSM:$sample_id\tLB:lib1" -o ${Name}.sam $refgenome "${pair1[i]}" "${pair2[i]}"
  samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -i ${Name}.sam | samtools view -S -b - > $Name.RG.NoDup.bam
  bam=$Name.RG.NoDup.bam
  # Extract the discordant paired-end alignments.
  #samtools view -b -F 1294 $bam > ${bam%.bam}.discordants.unsorted.bam
  # Extract the split-read alignments
  samtools view -h $bam \
	| /home/fs71400/grootcrego/software/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${bam%.bam}.splitters.unsorted.bam
  # Sort both alignments
  samtools sort -o ${bam%.bam}.discordants.sorted.bam ${bam%.bam}.discordants.unsorted.bam
  samtools sort -o ${bam%.bam}.splitters.sorted.bam ${bam%.bam}.splitters.unsorted.bam
done`

I ran lumpy-express, which is a fast wrapper of lumpy:
`~/software/lumpy-sv/bin/lumpyexpress -B $bam -S $splitters -D $discordants -o $output`

# Calling SV with Manta

First step is configuration and setup:
`cd $wd
for bam in $input ; do
	echo $bam
	cov=``echo $bam | cut -d'.' -f 2`
	`$manta --bam $bam --referenceFasta $ref --runDir /gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/calling_SV/manta/manta_$cov
	cd $wd
done`
