import sys
from Bio import SeqIO
import os.path

fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
# Get orthogroup ID, which is the name of the fastafile
og_id = str(sys.argv[1])
og_id = og_id.strip(".fasta")
outputfilename="checklist_orthologs_startstopcodons.txt"
# Since we iterate over a number of files, I only want to print the header once and then append additional lines when the file already exists
if os.path.isfile('checklist_orthologs_startstopcodons.txt'):
	output=open(outputfilename,'a')
else:
	output=open(outputfilename,'w')
	firstline = "gene_id\t\tog_id\tspecies\tstartcodon\tstopcodon\tlength\n"
	output.write(firstline)

# Define functions to check whether a startcodon or stopcodon is present at the start and end of the fasta sequence
def has_startcodon(seq):
	present = seq.startswith("ATG")
	return(present)

def has_stopcodon(seq):
	present = seq.endswith(("TAG", "TAA", "TGA"))
	return(present)

# compile infortmation for each fasta sequence in one line
for fasta in fasta_sequences:
	name = fasta.description
	name_split = name.split(" | ")
	gene_id = name_split[1]
	if gene_id.startswith("Tfas"):
		species = "fasciculata"
	if gene_id.startswith("Tlei"):
		species = "leiboldiana"
	sequence = str(fasta.seq)
	startcodon = has_startcodon(sequence)
	stopcodon = has_stopcodon(sequence)
	length = len(sequence)
	line_to_write =gene_id+"\t"+og_id+"\t"+species+"\t"+str(startcodon)+"\t"+str(stopcodon)+"\t"+str(length)+"\n"
	output.write(line_to_write)
