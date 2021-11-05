#!/usr/bin/env python

import sys
import re
from scipy.stats.distributions import chi2

#--------------
# Define input and output
null_file = open(sys.argv[1])
alt_file = open(sys.argv[2])
outputfilename="pairwise_dnds.codeml.allgenes.output"
output=open(outputfilename,'w')
#--------------

#--------------
# Define function to calculate likelihood ratio
def likelihood_ratio(llnull, llalt):
    return(2*(llalt-llnull))

non_decimal = re.compile(r'[^\d.]+')
#--------------

#--------------
# Main code: extract all relevant information and compile into one file
firstline="OrthoID"+"\t"+"dN/dS"+"\t"+"dN"+"\t"+"dS"+"\t"+"lnL_null"+"\t"+"lnL_alt"+"\t"+"Lratio"+"\t"+"pvalue\n"
output.write(firstline)

for line in null_file.readlines():
    if line.startswith("CODONML"):
        pre_ortho = line.split("/")[-1]
        ortho = pre_ortho[0:9]
    if line.startswith("lnL ="):
        line=line.replace('\n','')
        loglikelihood_null= float(line.split("=")[1])

for line in alt_file.readlines():
    if line.startswith("lnL ="):
        line=line.replace('\n','')
        loglikelihood_alt= float(line.split("=")[1])
    elif line.startswith("t= "):
        line=line.replace('\n','')
        line=line.replace(" ", "")
        dnds = non_decimal.sub('', line.split("=")[4])
        dn = non_decimal.sub('', line.split("=")[5])
        ds = non_decimal.sub('', line.split("=")[6])

LR = likelihood_ratio(loglikelihood_null,loglikelihood_alt)
p = chi2.sf(LR, 1) # one degree of freedom
p = '%.5f' % p

line_to_write = ortho + "\t" + dnds + "\t" + dn + "\t" + ds + "\t" + str(loglikelihood_null) + "\t" + str(loglikelihood_alt) + "\t" + str(LR) + "\t" + str(p)
print(line_to_write)
output.write(line_to_write)
