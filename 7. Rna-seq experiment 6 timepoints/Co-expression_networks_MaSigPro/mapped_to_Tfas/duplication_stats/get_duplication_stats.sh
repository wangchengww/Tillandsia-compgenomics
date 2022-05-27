modules='../Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic.txt'
ortho_info='../orthogroup-info_Significant_Genes-Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.size-corrected.noplastid-mito_EXONIC.txt'
ortho_counts='../orthogroup-counts_Significant_Genes-Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.size-corrected.noplastid-mito_EXONIC.txt'

echo "Module total one_one larger_Tfas larger_Tlei equal uniq_Tfas uniq_Tlei" > duplication_stats.per-orthogroup_EXONIC.txt
for i in $modules; do
  grep -f $i $ortho_info > list.txt
  cut -f 7 list.txt > list.ID.txt
  grep -f list.ID.txt $ortho_counts  > og_list.txt
  total=`wc -l og_list.txt | cut -f 1 -d' '`
  oneOne=`awk 'BEGIN { FS = "\t" } ; ($4 == $3) && ($4 == 1) {print $0}' og_list.txt | wc -l | cut -f 1 -d' '`
  larger_Tfas=`awk 'BEGIN { FS = "\t" } ; ($4 < $3) && !($4 == 0) {print $0}' og_list.txt | wc -l | cut -f 1 -d' '`
  larger_Tlei=`awk 'BEGIN { FS = "\t" } ; ($4 > $3) && !($3 == 0) {print $0}' og_list.txt | wc -l | cut -f 1 -d' '`
  equal=`awk 'BEGIN { FS = "\t" } ; ($3 == $4) && !($3 == 0) && !($3 == 1) {print $0}' og_list.txt | wc -l | cut -f 1 -d' '`
  uniq_Tfas=`awk 'BEGIN { FS = "\t" } ; ($4 == 0) && !($3 == 0) {print $0}' og_list.txt | wc -l | cut -f 1 -d' '`
  uniq_Tlei=`awk 'BEGIN { FS = "\t" } ; ($3 == 0) && !($4 == 0) {print $0}' og_list.txt | wc -l| cut -f 1 -d' ' `
  oneOne_p=`echo "scale=2 ; $oneOne / $total" | bc`
  larger_Tfas_p=`echo "scale=2 ; $larger_Tfas / $total" | bc`
  larger_Tlei_p=`echo "scale=2 ; $larger_Tlei / $total" | bc`
  equal_p=`echo "scale=2 ; $equal / $total" | bc`
  uniq_Tfas_p=`echo "scale=2 ; $uniq_Tfas / $total" | bc`
  uniq_Tlei_p=`echo "scale=2 ; $uniq_Tlei / $total" | bc`
  file=$i
  n=${file%.*}
  n=${n##*-}
  echo $n
  echo $n $total $oneOne_p $larger_Tfas_p $larger_Tlei_p $equal_p $uniq_Tfas_p $uniq_Tlei_p >> duplication_stats.per-orthogroup_EXONIC.txt
  rm list.txt list.ID.txt og_list.txt
done
