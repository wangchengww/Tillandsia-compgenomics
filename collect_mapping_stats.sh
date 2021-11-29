#!/bin/bash

file=/gpfs/data/fs71400/grootcrego/RERENCES_TILLANDSIA/RNA_experiment_6timepoints/2_mapped/subsample_mappingtest/*Log.final.out
echo "Name total uniq mismatch multi unmapped_short" > mapping_stats.txt
for i in $file
do
  NAME="$(basename $i Log.final.out)"
  echo $NAME
  total="$(grep "Number of input reads" $i | cut -f 2)"
  uniq="$(grep "Uniquely mapped reads %" $i | cut -f 2)"
  mismatch="$(grep "Mismatch rate per base, %" $i | cut -f 2)"
  multi="$(grep "% of reads mapped to multiple loci" $i | cut -f 2)"
  unmapped="$(grep "% of reads unmapped: too short" $i | cut -f 2)"
  echo $NAME $total $uniq $mismatch $multi $unmapped >> mapping_stats.txt
done
