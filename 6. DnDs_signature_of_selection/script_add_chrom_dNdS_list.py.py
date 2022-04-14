#!/usr/bin/env python
import sys
#--------------
# Define input and output
dnds_file = open(sys.argv[1])
orthogroup_file = open(sys.argv[2])
outputfilename=sys.argv[1].replace(".output",".with-chrom.output")
output=open(outputfilename,'w')
#--------------

#--------------
# Make dictionary of chromosome per orthogroup / species
og_dict = {}
for line in orthogroup_file.readlines():
    line = line.replace('\n','')
    splitted_line = line.split('\t')
    OG_ID = splitted_line[6]
    if splitted_line[0].startswith("Tfasc"):
        species = "Tfas"
    elif splitted_line[0].startswith("Tlei"):
        species = "Tlei"
    chrom = splitted_line[1]
	start = splitted_line[2]
	end = splitted_line[3]
    if OG_ID not in og_dict.keys():
        og_dict[OG_ID] = [species, chrom, start, end]
    else:
        og_dict[OG_ID].extend([species,chrom, start, end])
#--------------

#--------------
firstline="chrom_Tfas"+"\t"+"start"+"\t"+"end"+"\t"+"chrom_Tlei"+"\t"+"start"+"\t"+"end"+"\t"+"OrthoID"+"\t"+"dN/dS"+"\t"+"dN"+"\t"+"dS"+"\t"+"lnL_null"+"\t"+"lnL_alt"+"\t"+"Lratio"+"\t"+"pvalue\n"
output.write(firstline)

# Print out chromosome along with dnds
for line1 in dnds_file.readlines()[1:]:
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t') # split the line regarding the tabulations
    og = splitted_line1[0]
    chrom_tfas = og_dict[og][1]
	start_tfas = og_dict[og][2]
	end_tfas = og_dict[og][3]
    chrom_tlei=og_dict[og][5]
	start_tlei = og_dict[og][6]
	end_tlei = og_dict[og][7]
    new_line=chrom_tfas+"\t"+start_tfas+"\t"+end_tfas+"\t"+chrom_tlei+"\t"+start_tlei+"\t"+end_tlei+"\t"+line1+"\n"
    output.write(new_line)
