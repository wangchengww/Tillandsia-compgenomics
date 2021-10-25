#!/usr/bin/env python

import sys

vcf_in = open(sys.argv[1])  # auto-detect input format
coverage = int(sys.argv[2])
outputfilename=sys.argv[1].replace(".vcf",".filtered.0.4covsupport.PassPrecise.nohom.vcf")
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
		if format.startswith("1/1") or format.startswith("0/1"):
			if filter == "PASS" and precise == "PRECISE":
				if splitline[2].startswith("BND"):
					pe=splitinfo[6]
					splitPE=pe.split('=')
					PE = int(splitPE[1])
					sr = splitinfo[13]
					splitSR = sr.split('=')
					SR = int(splitSR[1])
					threshold = coverage*0.4
					if SR >= threshold and PE >= threshold:
						vcf_out.write(line)
				else:
					pe=splitinfo[4]
					splitPE=pe.split('=')
					PE = int(splitPE[1])
					sr = splitinfo[12]
					splitSR = sr.split('=')
					SR = int(splitSR[1])
					threshold = coverage*0.4
					if SR >= threshold and PE >= threshold:
						vcf_out.write(line)
