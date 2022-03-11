#!/usr/bin/env python

import sys
import re
from scipy.stats.distributions import chi2

#--------------
# Define input and output
null_file = open(sys.argv[1])
alt_file = open(sys.argv[2])
species=str(sys.argv[3])
outputfilename="branch-specific_dnds.{}.codeml.output".format(species)
output=open(outputfilename,'a+')
#--------------

#--------------
# Define function to calculate likelihood ratio
def likelihood_ratio(llnull, llalt):
    return(2*(llalt-llnull))

non_decimal = re.compile(r'[^\d.]+')
#--------------

#--------------
# Main code: extract all relevant information and compile into one file
if len(open(outputfilename, 'r').readlines(  )) == 0:
    firstline="OrthoID"+"\t"+"w_Tfas"+"\t"+"w_Tlei" + "\t"+"lnL_null"+"\t"+"lnL_alt"+"\t"+"Lratio"+"\t"+"pvalue\n"
    output.write(firstline)
for line in null_file.readlines():
    if line.startswith("CODONML"):
        pre_ortho = line.split("/")[-1]
        ortho = pre_ortho[0:9]
    if line.startswith("lnL"):
        line=line.replace('\n','')
        loglikelihood_null= float(line.split("  ")[4])
for line in alt_file:
	if line.startswith("lnL"):
		line=line.replace('\n','')
		loglikelihood_alt= float(line.split("  ")[4])
	elif line.startswith("w ratios as labels for TreeView"):
		line = next(alt_file)
		line=line.replace('\n','').replace('(', '').replace(')','').replace(' ','').replace('#',',')
		line=line.split(',')
		wtfas=line[1]
		wtlei=line[3]

LR = likelihood_ratio(loglikelihood_null,loglikelihood_alt)
p = chi2.sf(LR, 1) # one degree of freedom
p = '%.5f' % p

line_to_write = ortho + "\t" + wtfas + "\t" + wtlei +  "\t" + str(loglikelihood_null) + "\t" + str(loglikelihood_alt) + "\t" + str(LR) + "\t" + str(p) + "\n"
output=open(outputfilename,'a+')
output.write(line_to_write)
