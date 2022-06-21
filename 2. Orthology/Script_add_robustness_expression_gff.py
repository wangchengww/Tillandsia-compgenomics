#!/usr/bin/env python
# coding: utf-8

import sys
import re

gff = open(sys.argv[1])
genes = open(sys.argv[2])
ortho_info = open(sys.argv[3])
expression = pd.read_table(sys.argv[4])
expression_exon = pd.read_table(sys.argv[5])
directory_fasta = sys.argv[6]
outputfilename=sys.argv[1].replace(".gff",".with_robustness_orthology.gff")
output=open(outputfilename,'w')
outputfilename2="Robustness_checklist_ALLGenes.txt"
output2=open(outputfilename2,'w')

def has_startcodon(seq):
    present = seq.startswith("ATG")
    return(present)

def multiple_stopcodon(seq):
    x = len(seq)/3
    dividable = (x - int(x) == 0)
    if dividable == True:
        n = 3 # chunk length
        codons = [seq[i:i+n] for i in range(0, len(seq), n)]
        from collections import Counter
        counts = Counter(codons)
        stopcodon_counts = counts["TAG"]+counts["TAA"]+counts["TGA"]
        if stopcodon_counts == 1:
            multiple_stopcodons = "one_stopcodon"
        elif stopcodon_counts > 1:
            multiple_stopcodons = "multiple_stopcodons"
        elif stopcodon_counts == 0:
            multiple_stopcodons = "no_stopcodon"
    else:
        multiple_stopcodons = "not dividable by 3"
    return(multiple_stopcodons)

#1st step: find orthogroup and determine robustness
# Rule: a gene is ROBUS, when:
# - The gene is expressed across all exons
# OR
# - The gene has a start and stopcodon

first_line = "Gene_id\torthogroup\trobust\n"
output2.write(first_line)

gene_dict = {}
for line1 in genes.readlines():
    gene = line1.replace('\n','') # rm the return carriage
    print(gene)
	i = 0
	for line in ortho_info:
		if re.search(gene, line):
			i = 1
			og = line.split('\t')[6]
		if i == 0:
			og = "NO_ORTHOLOGY"
        # Extract info on length, start and stop codon from fasta sequence
    fasta_sequences = SeqIO.parse(open(directory_fasta_Tfas+gene+":cds.fst"),'fasta')
    for fasta in fasta_sequences:
		sequence = str(fasta.seq)
		startcodon = has_startcodon(sequence)
		multiple_stop = multiple_stopcodon(sequence)
		length = len(sequence)
        # Extract info on exon number from gff
        exon_count = 0
    gff_Tfas = open(sys.argv[1])
    for line in gff_Tfas:
            if re.search(gene, line):
                line = line.replace('\n','') # rm the return carriage
                splitted_line = line.split('\t')
                if splitted_line[2] == "exon":
                    exon_count=exon_count+1
        # Extract info on expression, the gene is only marked as unexpressed in a species when counts are 0 everywhere. Modification 17.06.2022: from now on, exons will be marked as not expressed when the average CPM across saples is < 0.001 (1000 counts in total)
        line = expression[expression['Geneid'].str.contains(gene)]
        if line.empty:
            expressed = "not_expressed"
        else:
            sum_Tfas = int(line.iloc[:, 1:37].sum(axis=1))
            if sum_Tfas > 0:
                exons = expression_exon[expression_exon['Geneid'].str.contains(gene)]
                expressed_exons = 0
                for row in exons.iterrows():
                    avg = int(line.iloc[:, 1:37].mean(axis=1))
                    if avg >= 0.001:
                        expressed_exons = expressed_exons + 1
                expressed = int(expressed_exons)/exon_count
                if expressed == 0:
                    expressed = "not_expressed"
            else:
                expressed = "not_expressed"
	if expressed == 1:
		robust = "ROBUST"
	elif startcodon == True and stopcodon == "one_stopcodon" && length >= 150:
		robust = "ROBUST"
	else:
		robust = "NOT_ROBUST"

    gene_list = [og, length, startcodon, multiple_stop, expressed, robust]
    gene_dict[gene] = gene_list
	to_write = gene+"\t"+og+"\t"+length+"\t"+startcodon+"\t"+multiple_stop+"\t"+expressed+"\t"+robrobust+"\n"
	output2.write(to_write)

#2nd step, iterate over each line of the gff and check if the ID is in the description dictionnary. If yes, print line + functional descriptions, if no just print the original line
for line2 in gff.readlines():
	line2 = line2.replace('\n','') # rm the return carriage
	splitted_line2 = line2.split('\t') # split the line regarding the tabulations
	annotation=splitted_line2[8]  # note that in python, it starts from 0, so the 9th field is [8]
	annotation2 = re.split('; | :', annotation)
	ID=annotation2[0].split('=')[1]
	if splitted_line2[2] == "gene":
		gene_key = any(key.startswith(ID) for key in gene_dict)
		stats = gene_dict[gene_key]
		line_to_write=line2+";"+stats[0]+";"+stats[4]+";"+stats[5]
		output.write(line_to_write)
	else:
		stats = gene_dict[ID]
		line_to_write=line2+";"+stats[0]+";"+stats[4]+";"+stats[5]
		output.write(line_to_write)
