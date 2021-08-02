#!/usr/bin/env python

import sys

per_gene = open(sys.argv[1])
counts = open(sys.argv[2])
outputfilename=sys.argv[1].replace(".txt",".table.txt")
output=open(outputfilename,'w')

og_dict = {}  # I initiate the dictionary
for line1 in counts.readlines()[1:]:
   line1 = line1.replace('\n','') # rm the return carriage
   splitted_line1 = line1.split('\t', 1) # split the line only at first occurring tab
   og = splitted_line1[0]
   counts_line1 = splitted_line1[1]
 og_dict[og] = counts_line1

for line2 in per_gene.readlines():
 line2 = line2.replace('\n','') # rm the return carriage
 splitted_line2 = line2.split('\t')
 og_id = splitted_line2[1]
 if (og_dict.has_key(og_id)):
	line_to_print = line2+"\t"+og_dict[og_id]+"\n"
	output.write(line_to_print)
