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

I decided to subsample samples A,B and D for both Tlei (D is individual used for the assembly) and A,B and F for Tfas (F is individual for the assembly), at times 01, 05, 13 and 17. This results in a total of 24 samples used in the test: 2 species x 3 samples x 4 timepoints.

I extracted random reads from each sample to 10 % of its original size using seqkit sample in the bash script `subsample.sh`. These were then mapped to each genome using STAR.
