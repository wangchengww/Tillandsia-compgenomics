#!/usr/bin/env python
import sys
#-------------------------
# Read input
gff_file = open("test") # For now optimized for .intact.gff3
total_genome_length = int(sys.argv[2])
scaffold_lengths=open(sys.argv[3])
# Name output
outputfilename_perscaff=sys.argv[1].replace(".gff",".per_scaffold.abundances.summary")
output_perscaff=open(outputfilename_perscaff,'w')

#--------------------------

#--------------------------
# Define functions
def perc(rel, tot):
	perc=(int(rel)/int(tot)*100)
	perc='%.2f' % perc
	return perc
#--------------------------

#-------------------------- WHOLE GENOME
# Initialize variables
total_genic_length=0

# Store scaffold lengths
lengths_dict = {}
for line in scaffold_lengths.readlines():
	line=line.replace("\n", "")
	scaffold = line.split("\t")[0]
	length = line.split("\t")[1]
	lengths_dict[scaffold] = length

# Read each line of the GFF, obtain length, type and subtype of TE, store in dictionnary per scaffold for per scaffold statistics meanwhile accumulate whole genome stats
per_scaffold_dict = {}
previous_scaffold="Scaffold_0"
for line in gff_file.readlines():
	if line.startswith("#"):
		continue
	else:
		if "Parent=" not in line:
			scaffold = line.split("\t")[0]
			if scaffold != previous_scaffold:
				previous_start=0
				previous_end=0
				previous_scaffold = scaffold
			start = int(line.split("\t")[3])
			end = int(line.split("\t")[4])
			if start < previous_end:
			# There is overlap
				if end < previous_end:
				# Nested TE
					length = 0
					previous_start = start
				else:
				# Partial overlap, 1 is added to length to account for GFF format
					start = previous_end
					length = int(end) - int(start) + 1
					previous_start = start
					previous_end = end
			else:
			# No overlap
				length = int(end) - int(start) + 1
				previous_start = start
				previous_end = end
			key = (scaffold, lengths_dict[scaffold])
			if key in per_scaffold_dict:
				per_scaffold_dict[key].append(length)
			else:
				key = (scaffold, lengths_dict[scaffold])
				per_scaffold_dict[key] = length
			#Tally the length per-type for whole genome stats
			total_genic_length = total_genic_length + length

# Perform whole genome stats
perc_total=perc(total_genic_length,total_genome_length)

#-------------------------- PER SCAFFOLD
# Re-initialize values
total_genic_length=0

# Print header in output
firstline="Scaffold\tlength\tproportion_coding\n"
output_perscaff.write(firstline)

# Obtain per scaffold statistics by iterating over dictionnary established in the previous step. Values are reset after each scaffold.
for key,value in per_scaffold_dict.items():
	scaffold=key[0]
	total_scaffold_length=key[1]
	for i in value:
		length = i
		total_genic_length = total_genic_length + length

	perc_total=perc(total_genic_length,total_scaffold_length)

	# Write per-scaffold stats to output
	output_perscaff.write(scaffold+"\t"+str(total_scaffold_length)+"\t"+str(perc_total)+"\n")
	# Reset all values
	total_genic_length=0
