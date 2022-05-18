#!/usr/bin/env python

from Bio import SeqIO
import re
import pandas as pd

outputfilename="Checklist_curated_orthologs_Tfas-Tlei.txt"
output=open(outputfilename,'w')

ortho_info = open('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt')
expression = pd.read_table('/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/RNA_experiment_6timepoints/3_counts/counts.Tfas_Tlei_6_timepoints.exons.edited.forR.sum.txt')
directory_fasta_Tlei = '/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tlei_assembly/assembly_26_scaffolds/fasta_seqs_curated_orthologs/'
directory_fasta_Tfas = '/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tfas_assembly/assembly_25_scaffolds/fasta_seq_curated_orthologs_CDS/'

def has_startcodon(seq):
	present = seq.startswith("ATG")
	return(present)

def has_stopcodon(seq):
	present = seq.endswith(("TAG", "TAA", "TGA"))
	return(present)

first_line='gene_id\torthogroup\tCDS_length\texon_count\tstartcodon\tstopcodon\texpressed_Tfas\texpressed_Tlei\n'
output.write(first_line)

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
        line = expression[expression['Geneid'].str.contains(gene)]
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
        expressed_Tfas = "NA"
        expressed_Tlei = "NA"
    # Compile information
    line_to_write=gene+'\t'+og+'\t'+str(length)+'\t'+str(exon_count)+'\t'+str(startcodon)+'\t'+str(stopcodon)+'\t'+str(expressed_Tfas)+'\t'+str(expressed_Tlei)+'\n'
    output.write(line_to_write)
