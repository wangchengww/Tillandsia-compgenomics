#!/usr/bin/env python

import sys


ogroups = open(sys.argv[1])
outputfilename=sys.argv[1].replace(".txt",".per_gene.txt")
output=open(outputfilename,'w')

for line in ogroups:
    line = line.replace('\r\n','')
    splitted_line = line.split('\t')
    OG_ID = splitted_line[1]
    genes_in_OG = splitted_line[4]+", "+splitted_line[5]
    genes_in_OG2 = genes_in_OG.split(", ")
    for gene in genes_in_OG2:
        if gene == '':
            continue
        else:
            line_to_print = gene+"\t"+OG_ID+"\n"
            output.write(line_to_print)
