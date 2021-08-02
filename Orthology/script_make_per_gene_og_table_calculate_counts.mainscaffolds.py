import sys
ogroups = open(sys.argv[1])
outputfilename=sys.argv[1].replace(".txt",".per_gene.txt")
output=open(outputfilename,'w')
og_dict = {}
for line in ogroups:
    line = line.replace('\r\n','')
    splitted_line = line.split('\t')
    OG_ID = splitted_line[1]
    genes_in_OG =splitted_line[3]+", "+splitted_line[4]+", "+splitted_line[5]
    genes_in_OG2 = genes_in_OG.split(", ")
    if OG_ID not in og_dict.keys():
        og_dict[OG_ID] = genes_in_OG
    else:
        og_dict[OG_ID] = og_dict[OG_ID]+", "+genes_in_OG

for key in og_dict:
    OG_ID = key
    genes_in_OG2 = og_dict[key]
    genes_in_OG2 = genes_in_OG2.split(", ")
    Acom = 0
    Tfas = 0
    Tlei = 0
    for gene in genes_in_OG2:
        if gene.startswith("Aco"):
            Acom = Acom + 1
        if gene.startswith("Tfasc"):
            Tfas = Tfas + 1
        if gene.startswith("Tlei"):
            Tlei = Tlei + 1
    for gene in genes_in_OG2:
        if gene == '':
            continue
        else:
            Acom = str(Acom)
            Tfas = str(Tfas)
            Tlei = str(Tlei)
            line_to_print = gene+"\t"+OG_ID+"\t"+Acom+"\t"+Tfas+"\t"+Tlei+"\n"
            print(line_to_print)
