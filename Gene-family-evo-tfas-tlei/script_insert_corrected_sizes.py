#!/usr/bin/env python

import sys

orthology_table = open(sys.argv[1])
corrected_sizes_tfas = open(sys.argv[2])
corrected_sizes_tlei = open(sys.argv[3])
outputfilename=sys.argv[1].replace(".txt",".size_corrections.txt")
output=open(outputfilename,'w')

og_dict_tfas = {}  # I initiate the dictionary
for line1 in corrected_sizes_tfas.readlines()[1:]:
	line1 = line1.replace('\n','') # rm the return carriage
	splitted_line1 = line1.split('\t') # split the line only at first occurring tab
	og = splitted_line1[0]
	counts_line1 = splitted_line1[2] # new family sizes
	og_dict_tfas[og] = counts_line1

og_dict_tlei = {}  # I initiate the dictionary
for line2 in corrected_sizes_tlei.readlines()[1:]:
	line2 = line2.replace('\n','') # rm the return carriage
	splitted_line2 = line2.split('\t') # split the line only at first occurring tab
	og = splitted_line2[0]
	counts_line2 = splitted_line2[2] # new family sizes
	og_dict_tlei[og] = counts_line2

for line3 in orthology_table.readlines():
	line3 = line3.replace('\n','') # rm the return carriage
	splitted_line3 = line3.split('\t')
	og_id = splitted_line3[6]
	if (og_id in og_dict_tfas) and (og_id in og_dict_tlei):
		line_to_print = ('\t'.join(str(x) for x in splitted_line3[0:8]))+"\t"+og_dict_tfas[og_id]+"\t"+ og_dict_tlei[og_id]+"\n"
	elif (og_id in og_dict_tfas):
		line_to_print=('\t'.join(str(x) for x in splitted_line3[0:8]))+"\t"+og_dict_tfas[og_id]+"\t"+ splitted_line3[9]+"\n"
	elif (og_id in og_dict_tlei):
		line_to_print=('\t'.join(str(x) for x in splitted_line3[0:9]))+"\t"+og_dict_tlei[og_id]+"\n"
	else:
		line_to_print=line3+"\n"
	output.write(line_to_print)
