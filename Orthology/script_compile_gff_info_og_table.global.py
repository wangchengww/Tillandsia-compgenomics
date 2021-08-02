#!/usr/bin/env python

import sys

orthogroups = open(sys.argv[1])
gff = open(sys.argv[2])
outputfilename=sys.argv[3]
output=open(outputfilename,'w')

# 1st, make a dictionnary of the gff with the information we want to transfer to the final table (scaffold, end and start position)

gff_dict = {}
for line1 in gff.readlines():
  line1 = line1.replace('\n','') # rm the return carriage
  splitted_line1 = line1.split('\t') # split the line regarding the tabulations
  if splitted_line1[2] == "mRNA":
	  pre_ID = splitted_line1[8]
	  pre_ID = pre_ID.split(";")
	  ID = pre_ID[0]
	  ID = ID.replace('ID=', '')
	  info = splitted_line1[0]+"\t"+splitted_line1[3]+"\t"+splitted_line1[4]
	  gff_dict[ID]=info
  else:
	  continue

# 2nd, iterate over the orthogroup list and print each line together with the gff information into a new table

for line2 in orthogroups.readlines():
  line2 = line2.replace('\n','') # rm the return carriage
  splitted_line2 = line2.split('\t') # split the line regarding the tabulations
  ID_OG = splitted_line2[0]
  if gff_dict.has_key(ID_OG):
	  line_to_print = ID_OG+"\t"+gff_dict[ID_OG]+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\n"
	  output.write(line_to_print)
  else:
	  line_to_print = ID_OG+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\n"
	  output.write(line_to_print)
