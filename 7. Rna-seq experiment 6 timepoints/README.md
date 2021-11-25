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

The final data was then moved to the VSC4 for further analyses.
