#!/usr/bin/env python

import sys

vcf_in = open(sys.argv[1])  # auto-detect input format
outputfilename=sys.argv[1].replace(".vcf",".filtered.10support.PassPrecise.vcf")
vcf_out=open(outputfilename,'w')

for line in vcf_in.readlines():
	if line.startswith('#'):
		vcf_out.write(line)
	else:
		splitline = line.split('\t')
		info=splitline[7]
		splitinfo=info.split(';')
		format=splitline[9]
		filter=splitline[6]
		precise=splitinfo[0]
		read_support=splitinfo[14]
		splitRE=read_support.split('=')
		RE = int(splitRE[1])
		if filter != "UNRESOLVED" and precise != "IMPRECISE" and RE >= 10:
			vcf_out.write(line)
