#!/usr/bin/env python

import sys
import re
from scipy.stats.distributions import chi2

#--------------
# Define input and output
null_file = open(sys.argv[1])
alt_file = open(sys.argv[2])
species=str(sys.argv[3])
outputfilename="branch-site_dnds.{}.codeml.output".format(species)
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
    firstline="OrthoID"+"\t"+"foreground_w"+"\t"+"2a_proportion"+"\t"+"lnL_null"+"\t"+"lnL_alt"+"\t"+"Lratio"+"\t"+"pvalue\n"
    output.write(firstline)
for line in null_file.readlines():
    if line.startswith("CODONML"):
        pre_ortho = line.split("/")[-1]
        ortho = pre_ortho[0:9]
    if line.startswith("lnL"):
        line=line.replace('\n','')
        loglikelihood_null= float(line.split("  ")[4])
for line in alt_file.readlines():
    if line.startswith("lnL"):
        line=line.replace('\n','')
        loglikelihood_alt= float(line.split("  ")[4])
    elif line.startswith("proportion"):
        line=line.replace('\n','')
        prop=list(filter(None, line.split(" ")))[3]
    elif line.startswith("foreground w"):
        line=line.replace('\n','')
        dnds = list(filter(None, line.split(" ")))[4]

LR = likelihood_ratio(loglikelihood_null,loglikelihood_alt)
p = chi2.sf(LR, 1) # one degree of freedom
p = '%.5f' % p

line_to_write = ortho + "\t" + dnds + "\t" + prop +  "\t" + str(loglikelihood_null) + "\t" + str(loglikelihood_alt) + "\t" + str(LR) + "\t" + str(p) + "\n"
output=open(outputfilename,'a+')
output.write(line_to_write)
