#!/usr/bin/env python

import sys

orthogroups = open(sys.argv[1])
scaffold = open(sys.argv[2])
gff_Tlei = open(sys.argv[2])
gff_Tfas = open(sys.argv[3])
gff_Acom = open(sys.argv[4])
outputfilename=sys.argv[5]
output=open(outputfilename,'w')

# 2nd, make a dictionnary of the gff with the information we want to transfer to the final table (scaffold, end and start position)

gff_dict_Tlei = {}
for line1 in gff_Tlei.readlines():
   line1 = line1.replace('\n','') # rm the return carriage
   splitted_line1 = line1.split('\t') # split the line regarding the tabulations
   if splitted_line1[2] == "mRNA":
	   pre_ID = splitted_line1[8]
	   pre_ID = pre_ID.split(";")
	   ID = pre_ID[0]
	   ID = ID.replace('ID=', '')
	   description_line = pre_ID[7]
	   GO = "NA"
	   for element in pre_ID:
		   if element.startswith("Ontology_id"):
			   ont = element
			   ont_split = ont.split('=')
			   GO = ont_split[1]
		   else:
			   continue
	   info = splitted_line1[0]+"\t"+splitted_line1[3]+"\t"+splitted_line1[4]+"\t"+GO+"\t"+description_line
	   gff_dict_Tlei[ID]=info


gff_dict_Tfas = {}
for line1 in gff_Tfas.readlines():
  line1 = line1.replace('\n','') # rm the return carriage
  splitted_line1 = line1.split('\t') # split the line regarding the tabulations
  if splitted_line1[2] == "mRNA":
	  pre_ID = splitted_line1[8]
	  pre_ID = pre_ID.split(";")
	  ID = pre_ID[0]
	  ID = ID.replace('ID=', '')
	  description_line = pre_ID[7]
	  GO = "NA"
	  for element in pre_ID:
		  if element.startswith("Ontology_id"):
			  ont = element
			  ont_split = ont.split('=')
			  GO = ont_split[1]
	  info = splitted_line1[0]+"\t"+splitted_line1[3]+"\t"+splitted_line1[4]+"\t"+GO+"\t"+description_line
	  gff_dict_Tfas[ID]=info
  else:
	  continue

gff_dict_Acom = {}
for line1 in gff_Acom.readlines():
  line1 = line1.replace('\n','') # rm the return carriage
  splitted_line1 = line1.split('\t') # split the line regarding the tabulations
  ID = splitted_line1[0]
  description_line = splitted_line1[3]
  GO = "NA"
  for element in splitted_line1:
	  if element.startswith("GO"):
		  GO = element
  position_line = splitted_line1[1]
  position_line = position_line.replace(':', '-')
  position_line = position_line.split('-')
  info = position_line[0]+"\t"+position_line[1]+"\t"+position_line[2]+"\t"+GO+"\t"+description_line
  gff_dict_Acom[ID]=info


# 3d, iterate over the orthogroup table and add scaffold length to each line

for line2 in orthogroups.readlines():
  line2 = line2.replace('\n','') # rm the return carriage
  splitted_line2 = line2.split('\t') # split the line regarding the tabulations
  ID_OG = splitted_line2[0]
  if gff_dict_Tlei.has_key(ID_OG):
	  line_to_print = ID_OG+"\t"+gff_dict_Tlei[ID_OG]+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\t"+splitted_line2[4]+"\n"
	  output.write(line_to_print)
  elif gff_dict_Tfas.has_key(ID_OG):
	  line_to_print = ID_OG+"\t"+gff_dict_Tfas[ID_OG]+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\t"+splitted_line2[4]+"\n"
	  output.write(line_to_print)
  elif gff_dict_Acom.has_key(ID_OG):
	  line_to_print = ID_OG+"\t"+gff_dict_Acom[ID_OG]+"\t"+splitted_line2[1]+"\t"+splitted_line2[2]+"\t"+splitted_line2[3]+"\t"+splitted_line2[4]+"\n"
	  output.write(line_to_print)
