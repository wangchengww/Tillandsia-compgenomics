#!/usr/bin/env python
import sys
#-------------------------
# Read input
gff_file = open(sys.argv[1]) # For now optimized for .intact.gff3
total_genome_length = int(sys.argv[2])
scaffold_lengths=open(sys.argv[3])
# Name output
outputfilename_full=sys.argv[1].replace(".gff3",".whole_genome.abundances.summary")
output_full=open(outputfilename_full,'w')
outputfilename_perscaff=sys.argv[1].replace(".gff3",".per_scaffold.abundances.summary")
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
total_repeat_length=0
cacta_length=0
mutator_length=0
pif_length=0
mariner_length=0
hat_length=0
helitron_length=0
copia_length=0
gypsy_length=0
unknown_length=0
mite_length=0
dna_length=0
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
			info= line.split("\t")[8]
			if line.split("\t")[2] == "repeat_region":
				type= (info.split(";")[2]).split("=")[1]
				subtype="none"
			else:
				type=line.split("\t")[2]
				subtype1=(info.split(";")[2]).split("=")[1]
				subtype=subtype1.split("/")[0]
				if type == "helitron":
					subtype = "none"
			list=[type,subtype,length]
			key = (scaffold, lengths_dict[scaffold])
			if key in per_scaffold_dict:
				per_scaffold_dict[key].append(list)
			else:
				key = (scaffold, lengths_dict[scaffold])
				per_scaffold_dict[key] = [list]
			#Tally the length per-type for whole genome stats
			total_repeat_length = total_repeat_length + length
			if type == "CACTA_TIR_transposon":
				cacta_length=cacta_length+length
			if type == "Mutator_TIR_transposon":
				mutator_length=mutator_length+length
			if type == "PIF_Harbinger_TIR_transposon":
				pif_length=pif_length+length
			if type == "Tc1_Mariner_TIR_transposon":
				mariner_length=mariner_length+length
			if type == "hAT_TIR_transposon":
				hat_length=hat_length+length
			if type == "helitron":
				helitron_length=helitron_length+length
			if type == "LTR/Copia":
				copia_length=copia_length+length
			if type == "Copia_LTR_retrotransposon":
				copia_length=copia_length+length
			if type == "LTR/Gypsy":
				gypsy_length=gypsy_length+length
			if type == "Gypsy_LTR_retrotransposon":
				gypsy_length=gypsy_length+length
			if type == "LTR/unknown":
				unknown_length=unknown_length+length
			if type == "LTR_retrotransposon":
				unknown_length=unknown_length+length
			if subtype == "MITE":
				mite_length=mite_length+length
			if subtype == "DNA":
				dna_length=dna_length+length


# Perform whole genome stats
perc_total=perc(total_repeat_length,total_genome_length)
perc_cacta=perc(cacta_length,total_genome_length)
perc_mutator=perc(mutator_length,total_genome_length)
perc_pif=perc(pif_length,total_genome_length)
perc_mariner=perc(mariner_length,total_genome_length)
perc_hat=perc(hat_length,total_genome_length)
perc_helitron=perc(helitron_length,total_genome_length)
perc_copia=perc(copia_length,total_genome_length)
perc_gypsy=perc(gypsy_length,total_genome_length)
perc_unknown=perc(unknown_length,total_genome_length)
perc_mite=perc(mite_length,total_genome_length)
perc_dna=perc(dna_length,total_genome_length)

# Write whole genome stats
output_full.write("Full Genome Annotation Report\n")
output_full.write("Total number of sequences: {}\n".format(len(per_scaffold_dict)))
output_full.write("Total number of bases: {}\n\n".format(total_genome_length))
output_full.write("TE type\t\tnumber of bp\tpercentage of bp\n")
output_full.write("Total\t\t{}\t{}\n".format(total_repeat_length,perc_total))
output_full.write("--- LTR ---\t--\t--\n")
output_full.write("Copia\t\t{}\t{}\n".format(copia_length,perc_copia))
output_full.write("Gypsy\t\t{}\t{}\n".format(gypsy_length, perc_gypsy))
output_full.write("Unknown\t\t{}\t{}\n".format(unknown_length, perc_unknown))
output_full.write("# # # Total:\t{}\t{}\n".format(int(gypsy_length)+int(copia_length)+int(unknown_length),float(perc_gypsy)+float(perc_copia)+float(perc_unknown)))
output_full.write("--- TIR ---\t--\t--\n")
output_full.write("CACTA\t\t{}\t\t{}\n".format(cacta_length,perc_cacta))
output_full.write("Mutator\t\t{}\t{}\n".format(mutator_length, perc_mutator))
output_full.write("pif_Harbinger\t{}\t\t{}\n".format(pif_length, perc_pif))
output_full.write("Tc1_Marinert\t{}\t\t{}\n".format(mariner_length, perc_mariner))
output_full.write("hAT\t\t{}\t\t{}\n".format(hat_length, perc_hat))
output_full.write("# # # Total:\t{}\t{}\n".format(int(cacta_length)+int(mutator_length)+int(pif_length)+int(mariner_length)+int(hat_length),float(perc_cacta)+float(perc_mutator)+float(perc_pif)+float(perc_mariner)+float(perc_hat)))
output_full.write("# # # DNA_transposons\t{}\t{}\n".format(dna_length,perc_dna))
output_full.write("# # # MITE\t\t{}\t\t{}\n".format(mite_length,perc_mite))
output_full.write("--- Helitron ---\t--\t--\n")
output_full.write("Helitron\t{}\t{}\n".format(helitron_length, perc_helitron))
#--------------------------

#-------------------------- PER SCAFFOLD
# Re-initialize values
total_repeat_length=0
cacta_length=0
mutator_length=0
pif_length=0
mariner_length=0
hat_length=0
helitron_length=0
copia_length=0
gypsy_length=0
unknown_length=0
mite_length=0
dna_length=0

# Print header in output
firstline="Scaffold\tlength\tRepetitive_content\tLTR/Copia\tLTR/Gypsy\tLTR/unknown\tTIR/CACTA\tTIR/Mutator\tTIR/PIF_Harbinger\tTIR/Tc1_Mariner\tTIR/hAT\tTIR/DNA_transposon\tTIR/MITE\tHelitron\n"
output_perscaff.write(firstline)

# Obtain per scaffold statistics by iterating over dictionnary established in the previous step. Values are reset after each scaffold.
for key,value in per_scaffold_dict.items():
	scaffold=key[0]
	total_scaffold_length=key[1]
	for list in value:
		type = list[0]
		subtype = list[1]
		length = list[2]
		total_repeat_length = total_repeat_length + length
		if type == "CACTA_TIR_transposon":
			cacta_length=cacta_length+length
		if type == "Mutator_TIR_transposon":
			mutator_length=mutator_length+length
		if type == "PIF_Harbinger_TIR_transposon":
			pif_length=pif_length+length
		if type == "Tc1_Mariner_TIR_transposon":
			mariner_length=mariner_length+length
		if type == "hAT_TIR_transposon":
			hat_length=hat_length+length
		if type == "helitron":
			helitron_length=helitron_length+length
		if type == "LTR/Copia":
			copia_length=copia_length+length
		if type == "Copia_LTR_retrotransposon":
			copia_length=copia_length+length
		if type == "LTR/Gypsy":
			gypsy_length=gypsy_length+length
		if type == "Gypsy_LTR_retrotransposon":
			gypsy_length=gypsy_length+length
		if type == "LTR/unknown":
			unknown_length=unknown_length+length
		if type == "LTR_retrotransposon":
			unknown_length=unknown_length+length
		if subtype == "MITE":
			mite_length=mite_length+length
		if subtype == "DNA":
			dna_length=dna_length+length
	perc_total=perc(total_repeat_length,total_scaffold_length)
	perc_cacta=perc(cacta_length,total_scaffold_length)
	perc_mutator=perc(mutator_length,total_scaffold_length)
	perc_pif=perc(pif_length,total_scaffold_length)
	perc_mariner=perc(mariner_length,total_scaffold_length)
	perc_hat=perc(hat_length,total_scaffold_length)
	perc_helitron=perc(helitron_length,total_scaffold_length)
	perc_copia=perc(copia_length,total_scaffold_length)
	perc_gypsy=perc(gypsy_length,total_scaffold_length)
	perc_unknown=perc(unknown_length,total_scaffold_length)
	perc_mite=perc(mite_length,total_scaffold_length)
	perc_dna=perc(dna_length,total_scaffold_length)
	# Write per-scaffold stats to output
	output_perscaff.write(scaffold+"\t"+str(total_scaffold_length)+"\t"+str(perc_total)+"\t"+str(perc_copia)+"\t"+str(perc_gypsy)+"\t"+str(perc_unknown)+"\t"+str(perc_cacta)+"\t"+str(perc_mutator)+"\t"+str(perc_pif)+"\t"+str(perc_mariner)+"\t"+str(perc_hat)+"\t"+str(perc_dna)+"\t"+str(perc_mite)+"\t"+str(perc_helitron)+"\n")
	# Reset all values
	total_repeat_length=0
	cacta_length=0
	mutator_length=0
	pif_length=0
	mariner_length=0
	hat_length=0
	helitron_length=0
	copia_length=0
	gypsy_length=0
	unknown_length=0
	mite_length=0
	dna_length=0
