import sys
import statistics

coverage = open(sys.argv[1]) # coverage file
outputfilename="Tfas_pergene_mediancov_and_orthoinfo.txt"
output=open(outputfilename,'w')

coverage = open("test_cov.txt")

# Calculate median coverage over each feature, and store in dictionary
median_dict = {}
prev_gene_id = ""
prev_other_gene = ""
gene_id = ""
other_gene = ""
covlist = []
overlap_covlist = []
i = 0
for line1 in coverage.readlines():
    i = i+1
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t') # split the line only at first occurring tab
    gene_id = splitted_line1[3]
    gene_id_spl = gene_id.split(',')
    if len(gene_id_spl) == 1:
        gene_id = splitted_line1[3]
        chrom = splitted_line1[0]
        pos = int(splitted_line1[1])
        cov = int(splitted_line1[2])
        if prev_gene_id == "":
            prev_gene_id = gene_id
            covlist.append(cov)
        elif gene_id == prev_gene_id:
            covlist.append(cov)
        elif gene_id != prev_gene_id:
            median_cov = statistics.median(covlist)
            mean_cov = statistics.mean(covlist)
            median_dict[prev_gene_id] = [str(median_cov),str(mean_cov),str(len(covlist))]
            print(prev_gene_id+"\t"+str(len(covlist)))
            covlist = []
            if gene_id == prev_other_gene:
                covlist = overlap_covlist
                overlap_covlist = []
            covlist.append(cov)
            prev_gene_id = gene_id
    elif len(gene_id_spl) > 1:
        gene_id = gene_id_spl[0]
        other_gene = gene_id_spl[1]
        chrom = splitted_line1[0]
        pos = int(splitted_line1[1])
        cov = int(splitted_line1[2])
        if prev_other_gene == "":
            prev_other_gene = other_gene
            overlap_covlist.append(cov)
        elif prev_other_gene == other_gene:
            overlap_covlist.append(cov)
        elif prev_other_gene != other_gene:
            if prev_other_gene in median_dict:
                prev_other_gene = other_gene
                overlap_covlist.append(cov)
            elif len(overlap_covlist) != 0:
                print(prev_other_gene+"\t"+str(len(overlap_covlist)))
                median_cov = statistics.median(overlap_covlist)
                mean_cov = statistics.mean(overlap_covlist)
                median_dict[prev_other_gene] = [str(median_cov),str(mean_cov),str(len(overlap_covlist))]
                overlap_covlist = []
                overlap_covlist.append(cov)
                prev_other_gene = other_gene
            else:
                overlap_covlist.append(cov)
                prev_other_gene = other_gene
        if prev_gene_id == "":
            prev_gene_id = gene_id
            covlist.append(cov)
        elif gene_id == prev_gene_id:
            covlist.append(cov)
        elif gene_id != prev_gene_id:
            median_cov = statistics.median(covlist)
            mean_cov = statistics.mean(covlist)
            median_dict[prev_gene_id] = [str(median_cov),str(mean_cov),str(len(covlist))]
            print(prev_gene_id+"\t"+str(len(covlist)))
            covlist = []
            covlist.append(cov)
            prev_gene_id = gene_id

median_dict[prev_gene_id] = [str(median_cov),str(mean_cov),str(len(covlist))]
print(prev_gene_id+"\t"+str(len(covlist)))
if prev_other_gene not in median_dict:
    median_dict[prev_other_gene] = [str(median_cov),str(mean_cov),str(len(overlap_covlist))]
    print(prev_other_gene+"\t"+str(len(overlap_covlist)))

firstline="gene_id\tmedian_cov\tmean_cov\tog_id\tAcom_count\tTfas_count\tTlei_count\n"
output.write(firstline)

ortho_info = open(sys.argv[2])

for line2 in ortho_info.readlines():
    line2 = line2.replace('\n','') # rm the return carriage
    splitted_line2 = line2.split('\t') # split the line only at first occurring tab
    gene_id = splitted_line2[0]
    if gene_id in median_dict:
        line_to_print=gene_id+"\t"+median_dict[gene_id][0]+"\t"+median_dict[gene_id][1]+"\t"+splitted_line2[6]+"\t"+splitted_line2[7]+"\t"+splitted_line2[8]+"\t"+splitted_line2[9]+"\n"
        output.write(line_to_print)
