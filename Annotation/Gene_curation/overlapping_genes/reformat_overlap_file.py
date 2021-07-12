#!/usr/bin/env python

import sys

overlap = open(sys.argv[1]) # bedfile
orthogroups = open(sys.argv[2])

outputfilename=sys.argv[1].replace(".txt",".table.txt")
output=open(outputfilename,'w')

# Store chr, stat and end position of each gene in dictionnary
gene_dict = {}  # I initiate the dictionary
for line1 in orthogroups.readlines():
	line1 = line1.replace('\n','') # rm the return carriage
	splitted_line1 = line1.split('\t')
	gene_ID0 = splitted_line1[0]
	OG1 = splitted_line1[6]
	gene_dict[gene_ID0] = OG1

first_line="gene_ID1\tOG_id1\tchr\tstart\tend\tlength\tstr\tdescription1\teAED1\tQI1\tgene_ID2\tOG_id2\tstart\tend\tstr\tdescription2\teAED2\tQI2\toverlap\n"
output.write(first_line)

for line1 in overlap.readlines():
	line1 = line1.replace('\n','') # rm the return carriage
	splitted_line1 = line1.split('\t') # split the line only at first occurring tab
	des1=splitted_line1[8]
	des1_spl=des1.split(';')
	gene_ID1=des1_spl[0].split('=')[1]
	if des1_spl[7].startswith('Description='):
		description1=des1_spl[7].split('=')[1]
	else:
		description1=des1_spl[7]
	length1 = int(splitted_line1[4]) - int(splitted_line1[3])
	des2=splitted_line1[17]
	des2_spl=des2.split(';')
	gene_ID2=des2_spl[0].split('=')[1]
	if des2_spl[7].startswith('Description='):
		description2=des2_spl[7].split('=')[1]
	else:
		description2=des2_spl[7]
	length2 = int(splitted_line1[13]) - int(splitted_line1[12])
	overlap = splitted_line1[-1]

	line_to_print = gene_ID1+'\t'+gene_dict[gene_ID1]+'\t'+splitted_line1[0]+'\t'+splitted_line1[3]+'\t'+splitted_line1[4]+'\t'+str(length1)+'\t'+splitted_line1[6]+'\t'+description1+'\t'+des1_spl[6].split('=')[1]+'\t'+des1_spl[5].split('=')[1]+'\t'+gene_ID2+'\t'+gene_dict[gene_ID2]+'\t'+splitted_line1[12]+'\t'+splitted_line1[13]+'\t'+str(length2)+'\t'+splitted_line1[15]+'\t'+description2+'\t'+des2_spl[6].split('=')[1]+'\t'+des2_spl[5].split('=')[1]+overlap+"\n"
	output.write(line_to_print)
