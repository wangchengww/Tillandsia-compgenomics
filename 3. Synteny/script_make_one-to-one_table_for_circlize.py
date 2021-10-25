import sys

orthogroups = open(sys.argv[1])
outputfilename=sys.argv[1].replace(".txt",".circlize.txt")
output=open(outputfilename,'w')

og_dict = {}  # I initiate the dictionary
for line1 in orthogroups.readlines():
	line1 = line1.replace('\n','') # rm the return carriage
	splitted_line1 = line1.split('\t') # split the line only at first occurring tab
	gene_id = splitted_line1[0]
	og = splitted_line1[6]
	if og in og_dict:
		list = [og_dict[og],gene_id+"\t"+splitted_line1[1]+"\t"+splitted_line1[2]+"\t"+splitted_line1[3]]
		og_dict[og] = list
	else:
		og_dict[og] = gene_id+"\t"+splitted_line1[1]+"\t"+splitted_line1[2]+"\t"+splitted_line1[3]

for key in og_dict:
	line_to_print = key+"\t"+og_dict[key][0]+"\t"+og_dict[key][1]+"\n"
	output.write(line_to_print)
