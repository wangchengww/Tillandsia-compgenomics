#!/usr/bin/env python

import sys

gff_in = open(sys.argv[1])
outputfilename=sys.argv[2]
output=open(outputfilename,'w')

for line in vcf_in.readlines():
	splitline = line.split('\t')
	info=splitline[8]
	splitinfo=info.split(';')
    genes = [i for i in splitinfo if i.startswith('Gene=')]
