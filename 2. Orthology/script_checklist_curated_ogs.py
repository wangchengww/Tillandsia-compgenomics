#!/usr/bin/env python

from Bio import SeqIO
import re
import pandas as pd

outputfilename="checklist_curated_orthologs_Tfas-Tlei.txt"

ortho_info = open('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/test')
expression_Tfas = pd.read_table('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/mapped_to_Tfas/counts.Tfas_Tlei_6_timepoints.exons.edited.forR.sum.txt')
expression_Tlei = pd.read_table('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/mapped_to_Tlei/counts.Tfas_Tlei_6_timepoints.exons.ToTLEI.edited.forR.sum.txt')
directory_fasta_Tlei = '/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tlei_assembly/assembly_26_scaffolds/fasta_seqs_curated_orthologs/'
directory_fasta_Tfas = '/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tfas_assembly/assembly_25_scaffolds/fasta_seq_curated_orthologs_CDS/'

def has_startcodon(seq):
	present = seq.startswith("ATG")
	return(present)

def has_stopcodon(seq):
	present = seq.endswith(("TAG", "TAA", "TGA"))
	return(present)

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
            stopcodon = has_stopcodon(sequence)
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
        # Extract info on expression, the gene is only marked as unexpressed in a species when counts are 0 everywhere
        line = expression_Tfas[expression_Tfas['Geneid'].str.contains(gene)]
        if line.empty:
            expressed_Tfas = "No features"
            expressed_Tlei = "No features"
        else:
            sum_Tfas = int(line.iloc[:, 1:37].sum(axis=1))
            sum_Tlei = int(line.iloc[:, 36:72].sum(axis=1))
            if sum_Tfas > 0:
                expressed_Tfas = True
            else:
                expressed_Tfas = False
            if sum_Tlei > 0:
                expressed_Tlei = True
            else:
                expressed_Tlei = False
    # Repeated operation for Tlei genes
    if gene.startswith("Tlei"):
        fasta_sequences = SeqIO.parse(open(directory_fasta_Tlei+gene+":cds.fst"),'fasta')
        for fasta in fasta_sequences:
            sequence = str(fasta.seq)
            startcodon = has_startcodon(sequence)
            stopcodon = has_stopcodon(sequence)
            length = len(sequence)
            exon_count=0
        gff_Tlei = open('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tlei_assembly/assembly_26_scaffolds/Tillandsia_leiboldiana_v1.2.edited_allfeatures.26chrom.curated-orthologs.gff')
        for line in gff_Tlei:
            if re.search(gene, line):
                line = line.replace('\n','') # rm the return carriage
                splitted_line = line.split('\t')
                if splitted_line[2] == "exon":
                    exon_count=exon_count+1
            # For now no expression information on Tlei genes
        line = expression_Tlei[expression_Tlei['Geneid'].str.contains(gene)]
        if line.empty:
            expressed_Tfas = "No features"
            expressed_Tlei = "No features"
        else:
            sum_Tfas = int(line.iloc[:, 1:37].sum(axis=1))
            sum_Tlei = int(line.iloc[:, 36:72].sum(axis=1))
            if sum_Tfas > 0:
                expressed_Tfas = True
            else:
                expressed_Tfas = False
            if sum_Tlei > 0:
                expressed_Tlei = True
            else:
                expressed_Tlei = False

    gene_list = [gene, og, length, exon_count, startcodon, stopcodon, expressed_Tfas, expressed_Tlei]
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
    df[8] = (max_exon_Tfas - df[3]).where(df[0].str.startswith("Tfas"), other=max_exon_Tlei - df[3])
    df[9] = (df[2]/max_length_Tfas).where(df[0].str.startswith("Tfas"), other=df[2]/max_length_Tlei)
    print(df)
    df.to_csv(outputfilename, index=False, header = False,sep='\t', mode='a')
