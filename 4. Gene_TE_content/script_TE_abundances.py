#!/usr/bin/env python

gff_file=open("test")
scaffold_lengths=open("test_lengths")

total_repeat_length=0
total_genome_length=837577910

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

def perc(rel, tot):
	perc=(int(rel)/int(tot)*100)
	perc='%.2f' % perc
	return perc

lengths_dict = {}
for line in scaffold_lengths.readlines():
	line=line.replace("\n", "")
	scaffold = line.split("\t")[0]
	length = line.split("\t")[1]
	lengths_dict[scaffold] = length

# Read each line of the GFF, obtain length, type and subtype of TE, store in dictionnary per scaffold for per scaffold statistics
per_scaffold_dict = {}
for line in gff_file.readlines():
	if "Parent=" not in line:
		scaffold = line.split("\t")[0]
		start = line.split("\t")[3]
		end = line.split("\t")[4]
		length = int(end) - int(start)
		info= line.split("\t")[8]
		if line.split("\t")[2] == "repeat_region":
			type= (info.split(";")[2]).split("=")[1]
			subtype="none"
		else:
			type=line.split("\t")[2]
			subtype=(info.split(";")[2]).split("=")[1]
			subtype=subtype.split("/")
		list=[type,subtype,length]
		key = (scaffold, lengths_dict[scaffold])
		if key in per_scaffold_dict:
			per_scaffold_dict[key].append(list)
		else:
			key = (scaffold, lengths_dict[scaffold])
			per_scaffold_dict[key] = [list]
		# Tally the length per-type for whole genome stats
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
		if type == "LTR/Gypsy":
			gypsy_length=gypsy_length+length
		if type == "LTR/unknown":
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

print("Full Genome Annotation Report\n")
print("Total number of sequences: {}\n".format(len(per_scaffold_dict)))
print("Total number of bases: {}\n\n".format(total_genome_length))
print("TE type\tnumber of bp\tpercentage of bp\n")
print("Total\t{}\t{}\n".format(total_repeat_length,perc_total))
print("--- LTR ---\n")
print("Copia\t{}\t{}\n".format(copia_length,perc_copia))
print("Gypsy\t{}\t{}\n".format(gypsy_length, perc_gypsy))
print("Unknown\t{}\t{}\n".format(unknown_length, perc_unknown))
print("--- TIR ---\n")
print("CACTA\t{}\t{}\n".format(cacta_length,perc_cacta))
print("Mutator\t{}\t{}\n".format(mutator_length, perc_mutator))
print("pif_Harbinger\t{}\t{}\n".format(pif_length, perc_pif))
print("Tc1_Marinert{}\t{}\n".format(mariner_length, perc_mariner))
print("hAT\t{}\t{}\n".format(hat_length, perc_hat))
print("\tDNA transposons\t{}\t{}\n".format(dna_length,perc_dna))
print("\tMITE\t{}\t{}\n".format(mite_length,perc_mite))
print("--- Helitron ---\n")
print("Helitron\t{}\t{}".format(helitron_length, perc_helitron))

# Obtain per-scaffold stats
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
		if type == "LTR/Gypsy":
			gypsy_length=gypsy_length+length
		if type == "LTR/unknown":
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
	print(scaffold+"\t"+str(total_scaffold_length)+"\t"+str(total_repeat_length)+"\t"+str(perc_total)+"\t"+str(copia_length)+"\t"+str(perc_copia))
