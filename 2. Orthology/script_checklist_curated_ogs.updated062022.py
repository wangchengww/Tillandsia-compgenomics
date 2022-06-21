#!/usr/bin/env python

# This script calculates whether there are multiple stop codons and shows expression as the percentage of exons expressed

from Bio import SeqIO
import re
import pandas as pd

outputfilename="checklist_curated_orthologs_Tfas-Tlei.new062022-2.txt"

ortho_info = open('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt')
expression_Tfas = pd.read_table('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/mapped_to_Tfas/counts.Tfas_Tlei_6_timepoints.exons.edited.forR.sum.txt')
Tfas_expression_exon = pd.read_table("/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/mapped_to_Tfas/counts.Tfas_Tlei_6_timepoints.exons.toTFAS.normalized-cpm.EdgeR.txt")
Tlei_expression_exon = pd.read_table("/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/mapped_to_Tlei/counts.Tfas_Tlei_6_timepoints.exons.toTLEI.normalized-cpm.EdgeR.txt")
expression_Tlei = pd.read_table('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/mapped_to_Tlei/counts.Tfas_Tlei_6_timepoints.exons.ToTLEI.edited.forR.sum.txt')
directory_fasta_Tlei = '/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tlei_assembly/assembly_26_scaffolds/fasta_seqs_curated_orthologs/'
directory_fasta_Tfas = '/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tfas_assembly/assembly_25_scaffolds/fasta_seq_curated_orthologs_CDS/'

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

og_dict = {}
for line1 in ortho_info.readlines():
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t')
    gene = splitted_line1[0]
    print(gene)
    og = splitted_line1[6]
    # Skip pineapple genes
    if gene.startswith("Aco"):
        continue
    if gene.startswith("Tfas"):
        # Extract info on length, start and stop codon from fasta sequence
        fasta_sequences = SeqIO.parse(open(directory_fasta_Tfas+gene+":cds.fst"),'fasta')
        for fasta in fasta_sequences:
            sequence = str(fasta.seq)
            startcodon = has_startcodon(sequence)
            multiple_stop = multiple_stopcodon(sequence)
            length = len(sequence)
            # Extract info on exon number from gff
            exon_count = 0
        gff_Tfas = open('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tfas_assembly/assembly_25_scaffolds/Tillandsia_fasciculata_v1.2.edited_allfeatures.25chrom.curated_orthologs_only.gff')
        for line in gff_Tfas:
            if re.search(gene, line):
                line = line.replace('\n','') # rm the return carriage
                splitted_line = line.split('\t')
                if splitted_line[2] == "exon":
                    exon_count=exon_count+1
        # Extract info on expression, the gene is only marked as unexpressed in a species when counts are 0 everywhere. Modification 17.06.2022: from now on, exons will be marked as not expressed when the average CPM across saples is < 0.001 (1000 counts in total)
        line = expression_Tfas[expression_Tfas['Geneid'].str.contains(gene)]
        if line.empty:
            expressed_Tfas = "No features"
            expressed_Tlei = "No features"
        else:
            sum_Tfas = int(line.iloc[:, 1:37].sum(axis=1))
            sum_Tlei = int(line.iloc[:, 36:72].sum(axis=1))
            if sum_Tfas > 0:
                exons = Tfas_expression_exon[Tfas_expression_exon['Geneid'].str.contains(gene)]
                expressed_exons = 0
                for row in exons.iterrows():
                    avg = int(line.iloc[:, 1:37].mean(axis=1))
                    if avg >= 0.001:
                        expressed_exons = expressed_exons + 1
                expressed_Tfas = int(expressed_exons)/exon_count
                if expressed_Tfas == 0:
                    expressed_Tfas = "not_expressed"
            else:
                expressed_Tfas = "not_expressed"
            if sum_Tlei > 0:
                exons = Tfas_expression_exon[Tfas_expression_exon['Geneid'].str.contains(gene)]
                expressed_exons = 0
                for row in exons.iterrows():
                    avg = int(line.iloc[:, 36:72].mean(axis=1))
                    if avg >= 0.001:
                        expressed_exons = expressed_exons + 1
                expressed_Tlei = int(expressed_exons)/exon_count
                if expressed_Tlei == 0:
                    expressed_Tlei = "not_expressed"
            else:
                expressed_Tlei = "not_expressed"
    # Repeated operation for Tlei genes
    if gene.startswith("Tlei"):
        fasta_sequences = SeqIO.parse(open(directory_fasta_Tlei+gene+":cds.fst"),'fasta')
        for fasta in fasta_sequences:
            sequence = str(fasta.seq)
            startcodon = has_startcodon(sequence)
            multiple_stop = multiple_stopcodon(sequence)
            length = len(sequence)
            exon_count=0
        gff_Tlei = open('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tlei_assembly/assembly_26_scaffolds/Tillandsia_leiboldiana_v1.2.edited_allfeatures.26chrom.curated-orthologs.gff')
        for line in gff_Tlei:
            if re.search(gene, line):
                line = line.replace('\n','') # rm the return carriage
                splitted_line = line.split('\t')
                if splitted_line[2] == "exon":
                    exon_count=exon_count+1
        line = expression_Tlei[expression_Tlei['Geneid'].str.contains(gene)]
        if line.empty:
            expressed_Tfas = "No features"
            expressed_Tlei = "No features"
        else:
            sum_Tfas = int(line.iloc[:, 1:37].sum(axis=1))
            sum_Tlei = int(line.iloc[:, 36:72].sum(axis=1))
            if sum_Tfas > 0:
                exons = Tlei_expression_exon[Tlei_expression_exon['Geneid'].str.contains(gene)]
                expressed_exons = 0
                for row in exons.iterrows():
                    avg = int(line.iloc[:, 1:37].mean(axis=1))
                    if avg >= 0.001:
                        expressed_exons = expressed_exons + 1
                expressed_Tfas = int(expressed_exons)/exon_count
                if expressed_Tfas == 0:
                    expressed_Tfas = "not_expressed"
            else:
                expressed_Tfas = "not_expressed"
            if sum_Tlei > 0:
                exons = Tlei_expression_exon[Tlei_expression_exon['Geneid'].str.contains(gene)]
                expressed_exons = 0
                for row in exons.iterrows():
                    avg = int(line.iloc[:, 36:72].mean(axis=1))
                    if avg >= 0.001:
                        expressed_exons = expressed_exons + 1
                expressed_Tlei = int(expressed_exons)/exon_count
                if expressed_Tlei == 0:
                    expressed_Tlei = "not_expressed"
            else:
                expressed_Tlei = "not_expressed"

    gene_list = [gene, og, length, exon_count, startcodon, multiple_stop, expressed_Tfas, expressed_Tlei]
    if og in og_dict:
        og_dict[og].append(gene_list)
    else:
        lst = []
        lst.append(gene_list)
        og_dict[og] = lst

column_names = ["Gene_id", "orthogroup", "CDS_length", "exon_count", "startcodon", "stopcodon", "expressed_Tfas", "expressed_Tlei", "diff_largest_exon_count", "diff_longest_CDS"]
df0 = pd.DataFrame(columns = column_names)
df0.to_csv(outputfilename, index=False, sep='\t', mode='w')

for og, genes in og_dict.items():
    df = pd.DataFrame(genes)
    max_exon_Tfas = df.loc[df[0].str.startswith("Tfas"), 3].max()
    max_exon_Tlei = df.loc[df[0].str.startswith("Tlei"), 3].max()
    max_length_Tfas = df.loc[df[0].str.startswith("Tfas"), 2].max()
    max_length_Tlei = df.loc[df[0].str.startswith("Tlei"), 2].max()
    df[9] = (max_exon_Tfas - df[3]).where(df[0].str.startswith("Tfas"), other=max_exon_Tlei - df[3])
    df[10] = (df[2]/max_length_Tfas).where(df[0].str.startswith("Tfas"), other=df[2]/max_length_Tlei)
    df.to_csv(outputfilename, index=False, header = False,sep='\t', mode='a')
