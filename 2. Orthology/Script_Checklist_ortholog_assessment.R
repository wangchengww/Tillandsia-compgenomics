setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/2. Orthology/")

data <- read.delim("checklist_curated_orthologs_Tfas-Tlei.txt", header = T, sep = "\t")
# Count number of genes in T. fasciculata with the same number of exons and the same 
# length as the largest gene in the orthogroup
dim(data[data$diff_largest_exon_count == 0 & 
           startsWith(data$Gene_id, "Tfas"),]) # 19738 (75 %)
dim(data[data$diff_longest_CDS == 1 & 
           startsWith(data$Gene_id, "Tfas"),]) # 18416 (69.9 % %)
# Count number of genes in T. leiboldiana with the same number of exons as the largest
# gene in the orthogroup
dim(data[data$diff_largest_exon_count == 0 & 
           startsWith(data$Gene_id, "Tlei"),]) # 19786 (84 %)
dim(data[data$diff_longest_CDS == 1 & 
           startsWith(data$Gene_id, "Tlei"),]) # 18353 (77.8 % %)

library(dplyr)
tfas_data <- data[startsWith(multicopy$Gene_id, "Tfas"),] 
tlei_data <- data[startsWith(multicopy$Gene_id, "Tlei"),]
tfas_multicopy <-  tfas_data %>% group_by(orthogroup) %>% filter(n()>2) # 9382 genes
dim(tfas_multicopy[tfas_multicopy$diff_largest_exon_count == 0,]) # 4197 (44.7 %)
dim(tfas_multicopy[tfas_multicopy$diff_longest_CDS == 1,]) # 3085 (32.8 %)

tlei_multicopy <-  tlei_data %>% group_by(orthogroup) %>% filter(n()>2) # 6193 genes
dim(tlei_multicopy[tlei_multicopy$diff_largest_exon_count == 0,]) # 2715 (43.8 %)
dim(tlei_multicopy[tlei_multicopy$diff_longest_CDS == 1,]) # 1969 (31.8 %)
# Count number of genes in T. leiboldiana with the same number of exons as the largest
# gene in the orthogroup
dim(multicopy[multicopy$diff_largest_exon_count == 0 & 
           startsWith(multicopy$Gene_id, "Tlei"),]) # 5398 (84 %)
dim(multicopy[multicopy$diff_longest_CDS == 1 & 
           startsWith(multicopy$Gene_id, "Tlei"),]) # 4087 (77.8 % %)
