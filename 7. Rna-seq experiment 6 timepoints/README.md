# 6 timepoint RNA-seq experiment between CAM and C3 Tillandsia

This folder reports the steps undertaken to infer differential gene expression across 6 time points within a 24 hour period between a CAM Tillandsia (*Tillandsia fasciculata*) and a C3 plant (*Tillandsia leiboldiana*).

# Raw data processing

The sequencing data were delivered to us already demultiplexed. Since two rounds of sequencing was done, I had to collapse the two parts of data into one for each sample. I first modified sample names:

  input_dir=/media/clara/data2/RNA_seq_Tlei_Tfas_6timepoints/raw/HT3NKDSX2_1_R12435_20211120/demultiplexed
  output_dir=/media/clara/data2/RNA_seq_Tlei_Tfas_6timepoints/combined/lane2
  cat change_names.txt | while read line; do
    Name=`echo "$line"|awk '{print $1}'`
    Replace=`echo "$line"|awk '{print $2}'`
    echo $Name
    cd $input_dir/$Name
    cp ${Name}*R1*.fastq.gz $output_dir/${Replace}.R1.fastq.gz
    cp ${Name}*R2*.fastq.gz $output_dir/${Replace}.R2.fastq.gz
  done

Then collapsed both parts of the data into one file:

  lane1=/media/clara/data2/RNA_seq_Tlei_Tfas_6timepoints/combined/lane1
  lane2=/media/clara/data2/RNA_seq_Tlei_Tfas_6timepoints/combined/lane2
  output_dir=/media/clara/data2/RNA_seq_Tlei_Tfas_6timepoints/combined/full_data
  cat change_names.txt | while read line; do
    basename=`echo "$line"|awk '{print $2}'`
    cat $lane1/${basename}.R1.fastq.gz $lane2/${basename}.R1.fastq.gz > $output_dir/${basename}.R1.fastq.gz
    cat $lane1/${basename}.R2.fastq.gz $lane2/${basename}.R2.fastq.gz > $output_dir/${basename}.R2.fastq.gz
  done

The final data was then moved to the VSC4 for further analyses. Trimming was performed using AdapterRemoval v.2.3.1 with the bash script `trim.sh`. Post-trimming, samples contained 58 - 97 million pairs, with three outlier samples containing about 150 million pairs. Our data contains a fair amount of overrepresented sequences in pair 2, of which many are simply a sequence of G's - probably a by-product of the sequencing approach.

# Choosing a reference genome for mapping

Since most of our RNA analysis will be based on counts of reads per gene, it is important to choose a reference genome that doesn't create inherent biases in the counts. Since our references are from *T. leiboldiana* and *T. fasciculata*, I am devising a small test where I map a subset of our data to both references and also to *A. comosus* to infer possible mapping bias and choose the better reference genome for this analysis.

I decided to subsample samples A,C and D for both Tlei (D is individual used for the assembly) and A,C and F for Tfas (F is individual for the assembly), at times 01, 05, 13 and 17. This results in a total of 24 samples used in the test: 2 species x 3 samples x 4 timepoints.

I extracted random reads from each sample to 10 % of its original size using seqkit sample in the bash script `subsample.sh`. These were then mapped to each genome using STAR with `map.sh`. Mapping statistics were collected from each log file with the bash script `collect_mapping_stats.sh` and assessed with the r script `Assessment_mapping_bias.R`. The assessment showed that RNA data maps far better to the Tillandsia genomes than to A. comosus (around 90 % to around 50 %). Mapping bias seems much stronger when mapping to T. leiboldiana than to T. fasciculata. This seems to be rather due to unmapped reads than due to multimapper reads. Based on these results, I decided to map the full dataset to T. fasciculata.

I repeated the test using the masked versions of both Tillandsia assemblies to see if this would perhaps show different rates of mapping bias. However, the mapping bias results are very similar, meaning that bias differences between T. fasciculata and T. leiboldiana is not due to active TEs is the RNA-seq.

# Mapping full dataset

All 72 samples were then mapped to the *T. fasciculata* genome with the same script as above. General mapping stats were obtained with MultiQC and also with the script above. Further assessment can be found in `Assessment_mapping_bias.R`

# Obtaining counts per gene

Counts per gene were computed using featureCounts with the bash script `counts.sh`. The summary file was then transposed using datamash: `datamash transpose < counts.Tfas_Tlei_6_timepoints.txt.summary > transposed_file`
I also computed the number of genes per sample with 0 counts:

	for i in {7..79}; do
	 cut -f $i counts.Tfas_Tlei_6_timepoints.txt | grep -c -w 0
	done

Generally, about 72 % - 78 % of alignments were assinged to a gene. 15 % - 20 % could not be aligned due to multiple mapping, 5 % - 7 % due to being outside any gene. 25 % - 30 % of the genes did not have any alignments per sample.

The count data was inspected with Principal Component Analysis in the R sheet `Inspect_count_data.R`. The PCA showed that samples are primarily separated by species (32 % variance explained). When making PCA's only for one or the other species, PC's tend to separate by individual (and not by time point), though some individuals seem more similar in counts than others (See google presentation of Drive for the PCA plots). Therefore, we have to count with much between-individual variation in time-dependent expression in downstream analyses. I repeated these analyses based only on genes called as orthologs with T.lei and A.com (~ 26,000 genes) to see if individual variation was somehow dependent on low-quality gene models or TEs. However, the PCA was identical, meaning that the per-individual variation that I observe in the count data is not based on "bad genes".

# Differential Gene expression

Within species comparisons, Nearby time points:

1 am - 5 am
5 am - 9 am
9 am - 1 pm
1 pm - 5 pm
5 pm - 9 pm
9 pm - 1 am

For each species, we get a list of DE genes. This will give us an idea of gradual changes in expression across time. We could construct a table at which genes are "on/off" by noting their changes in expression. It may also be worth it to do DGE at larger time points, as changes could be too gradual to pick up over short times. So, we could also do DGE over 8 hours:

1 am - 9 am
5 am - 1 pm
9 am - 5 pm
1 pm - 9 pm
5 pm - 1 am

The overlap could then be compared between species, or the timing at which certain genes are "on/off".

# Co-expression analysis

Starts from the count data. Data needs to be normalized, which is apparently a different normalization approach from classic DGE. Possible packages are TMM, DEseq2
...
Co-expression modules can be called using WCGNA (most classical approach), but DICER could also be interesting: "By design, DICER is tailored to identify module pairs that correlate differently between sample groups, e.g. modules that form one large interconnected module in one group compared with several smaller modules in another (Figure 3D). DICER may be particularly useful for time series ex- periments in which co-expression changes are gradual, e.g. cell cycle series experiments, where modules are specific to a par- ticular phase and co-expressed in transitions between phases"
...
Since WGCNA and DICER identify slightly different types of module changes, both approaches could be used and combined
