#!/usr/bin/env python

import sys

orthology_table = open(sys.argv[1])
corrected_sizes = open(sys.argv[2])
outputfilename=sys.argv[1].replace(".txt",".size_corrections.txt")
output=open(outputfilename,'w')

og_dict = {}  # I initiate the dictionary
for line1 in corrected_sizes.readlines()[1:]:
	line1 = line1.replace('\n','') # rm the return carriage
	splitted_line1 = line1.split('\t') # split the line only at first occurring tab
	og = splitted_line1[0]
	counts_line1 = splitted_line1[2] # new family sizes
	og_dict[og] = counts_line1

for line2 in orthology_table.readlines():
	line2 = line2.replace('\n','') # rm the return carriage
	splitted_line2 = line2.split('\t')
	og_id = splitted_line2[6]
	if og_id in og_dict:
		line_to_print = ('\t'.join(str(x) for x in splitted_line2[0:8]))+"\t"+og_dict[og_id]+"\t"+ splitted_line2[9]+"\n"
	else:
		line_to_print=line2+"\n"
	output.write(line_to_print)
