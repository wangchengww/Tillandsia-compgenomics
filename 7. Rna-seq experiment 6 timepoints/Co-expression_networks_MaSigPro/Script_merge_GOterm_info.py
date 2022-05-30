#!/usr/bin/env python
import sys

go_tfas= open(sys.argv[1])
go_tlei= open(sys.argv[2])
outputfilename=sys.argv[3]
output=open(outputfilename,'w')

go_dict = {}
for line1 in go_tfas.readlines()[1:]:
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t')
    term = splitted_line1[1]
    go_dict[term] = list( splitted_line1[i] for i in [0,1,2,3,5] )

for line2 in go_tlei.readlines()[1:]:
    line2 = line2.replace('\n','') # rm the return carriage
    splitted_line2 = line2.split('\t')
    term = splitted_line2[1]
    if term in go_dict:
        go_dict[term].append(splitted_line2[3])
        go_dict[term].append(splitted_line2[5])
    else:
        go_dict[term] = list( splitted_line2[i] for i in [0,1,2,3,5] )
firstline='Category\tID\tTerm\tadj_pval_tfas\tadj_pval_tlei\tGenes\n'
output.write(firstline)
for go, info in go_dict.items():
    if (len(info) > 5):
        genes=info[4]+', '+info[6]
        pval_tfas=info[3]
        pval_tlei=info[5]
        line=info[0]+'\t'+info[1]+'\t'+info[2]+'\t'+pval_tfas+'\t'+pval_tlei+'\t'+genes+'\n'
    else:
        genes=info[4]
        if genes.startswith("Tfas"):
            pval_tfas=info[3]
            line=info[0]+'\t'+info[1]+'\t'+info[2]+'\t'+pval_tfas+'\tNA\t'+genes+'\n'
        else:
            pval_tlei=info[3]
            line=info[0]+'\t'+info[1]+'\t'+info[2]+'\tNA\t'+pval_tlei+'\t'+genes+'\n'
    output.write(line)
