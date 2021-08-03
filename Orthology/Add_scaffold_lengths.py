#!/usr/bin/env python

import sys

orthogroups = open(sys.argv[1])
scaffold = open(sys.argv[2])
outputfilename=sys.argv[1].replace(".txt",".per_scaffold.txt")
output=open(outputfilename,'w')

# 1st, make a dictionnary of the scaffold ID (key) and length (value)

scaffold_dict = {}
for line1 in scaffold.readlines():
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t') # split the line regarding the tabulations
    ID = splitted_line1[0]
    length = splitted_line1[1]
    scaffold_dict[ID] = length

# 2nd, iterate over the orthogroup table and add scaffold length to each line

for line2 in orthogroups.readlines():
    line2 = line2.replace('\n','') # rm the return carriage
    splitted_line2 = line2.split('\t') # split the line regarding the tabulations
    scaffold_OG = splitted_line2[1]
    if scaffold_dict.has_key(scaffold_OG):
        line_to_print = splitted_line2[0]+"\t"+splitted_line2[1]+"\t"+scaffold_dict[scaffold_OG]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\t"+splitted_line2[6]+"\t"+splitted_line2[7]+"\t"+splitted_line2[8]+"\t"+splitted_line2[9]+"\t"+splitted_line2[4]+"\t"+splitted_line2[5]+"\n"
        output.write(line_to_print)
