#!/usr/bin/env python

import sys

vcf_in = open(sys.argv[1])  # auto-detect input format
coverage = int(sys.argv[2])
outputfilename=sys.argv[1].replace(".vcf",".filtered.0.4covsupport.Precise.nohom.vcf")
vcf_out=open(outputfilename,'w')

for line in vcf_in.readlines():
	if line.startswith('#'):
		vcf_out.write(line)
	else:
		line_nosp = line.replace('\n','')
		splitline = line_nosp.split('\t')
		info=splitline[7]
		splitinfo=info.split(';')
		format=splitline[9]
		if format.startswith("1/1") or format.startswith("0/1"):
			if "IMPRECISE" not in splitinfo:
				pe=splitinfo[-2]
				splitPE=pe.split('=')
				PE = int(splitPE[1])
				sr = splitinfo[-1]
				splitSR = sr.split('=')
				SR = int(splitSR[1])
				threshold = coverage*0.4
				if SR >= threshold and PE >= threshold:
					vcf_out.write(line)
