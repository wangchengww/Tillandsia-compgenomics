import sys

regions = open(sys.argv[1]) # bedfile

# Store chr, stat and end position of each gene in dictionnary
gene_dict = {}  # I initiate the dictionary
for line1 in regions.readlines():
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t') # split the line only at first occurring tab
    gene_id = splitted_line1[3]
    chrom = splitted_line1[0]
    start = int(splitted_line1[1])
    end = int(splitted_line1[2])
    for i in range(start, end):
        gene_dict[(chrom, i)] = [chrom,gene_id]

coverage = open(sys.argv[2]) # coverage file
outputfilename=sys.argv[2].replace(".txt",".edited.txt")
output=open(outputfilename,'w')

# Calculate median average
for line2 in coverage.readlines():
    line2 = line2.replace('\n','') # rm the return carriage
    splitted_line2 = line2.split('\t') # split the line only at first occurring tab
    chrom1 = splitted_line2[0]
    pos = int(splitted_line2[1])
    if (chrom1, pos) in gene_dict:
        gene_id = gene_dict[(chrom1,pos)][1]
        line_to_print = line2+"\t"+gene_id+"\n"
        output.write(line_to_print)
