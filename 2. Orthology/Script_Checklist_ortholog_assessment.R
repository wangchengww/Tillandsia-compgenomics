setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/2. Orthology/")
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/2. Orthology/")

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
library(ggplot2)
library(grid)
library(gridExtra)
tfas_data <- data[startsWith(data$Gene_id, "Tfas"),]
dim(tfas_data[tfas_data$startcodon == "True" & tfas_data$stopcodon == "True",]) # 23,618 (89.7 % of all genes are complete)
tlei_data <- data[startsWith(data$Gene_id, "Tlei"),]
dim(tlei_data[tlei_data$startcodon == "True" & tlei_data$stopcodon == "True",]) # 14,924 (63,3 % of all genes are complete)
tlei_incomplete <- tlei_data[!(tlei_data$startcodon == "True" & tlei_data$stopcodon == "True"),]
write.table(tlei_incomplete, file = "Incomplete_genes_start-stopcodon_Tlei.txt", 
            quote = F, row.names = F, sep = "\t")
 # How many multicopy genes are complete in exon count and length in Tfas
tfas_multicopy <-  tfas_data %>% group_by(orthogroup) %>% filter(n()>1) # 12014 genes
length(unique(tfas_multicopy$orthogroup)) # 3331 orthogroups, meaning that there are 8683 secondary copies
dim(tfas_multicopy[tfas_multicopy$diff_largest_exon_count == 0,]) # 5427 (45 %) of which 2096 (24 %) are a secondary copy
dim(tfas_multicopy[tfas_multicopy$diff_longest_CDS == 1,]) # 4105 (34.2%) of which 774 are a secondary copy

# How many multicopy genes are complete in exon count and length in Tlei
tlei_multicopy <-  tlei_data %>% group_by(orthogroup) %>% filter(n()>1) # 8011 genes
length(unique(tlei_multicopy$orthogroup)) # 2311 orthogroups, meaning that 5700 genes are secondary copies
dim(tlei_multicopy[tlei_multicopy$diff_largest_exon_count == 0,]) # 4213 (52.6 %) of which 1902 (34 %) are a secondary copy
dim(tlei_multicopy[tlei_multicopy$diff_longest_CDS == 1,]) # 2780 (34.7 %) of which 469 are a secondary copy

# How many multicopy Tfas genes are expressed in both species, and how does this compare to completeness
dim(tfas_multicopy[tfas_multicopy$expressed_Tfas == "True",]) # 9994 of which 6663 are a secondary copy
dim(tfas_multicopy[tfas_multicopy$expressed_Tfas == "True" & tfas_multicopy$diff_largest_exon_count == 0,]) # 4923 (90 % of complete genes are expressed in Tfas) 
dim(tfas_multicopy[tfas_multicopy$expressed_Tlei == "True" & tfas_multicopy$diff_largest_exon_count == 0,]) # 4844 (89 % of complete genes are expressed in Tlei)

# How many multicopy Tlei genes are expressed in both species
dim(tlei_multicopy[tlei_multicopy$expressed_Tfas == "True",]) # 6194 of which 3883 are a secondary copy
dim(tlei_multicopy[tlei_multicopy$expressed_Tlei == "True",]) # 6443 of which 4112 are a secondary copy

# Make histogram of exon count and length differences per family
p1 <-ggplot(data = tfas_multicopy, aes(x = diff_longest_CDS)) +
  geom_histogram()
p2 <- ggplot(data = tlei_multicopy, aes(x = diff_longest_CDS)) +
  geom_histogram()
grid.arrange(p1, p2, nrow = 1)

p1 <-ggplot(data = tfas_multicopy, aes(x = diff_largest_exon_count)) +
  geom_histogram()
p2 <- ggplot(data = tlei_multicopy, aes(x = diff_largest_exon_count)) +
  geom_histogram()
grid.arrange(p1, p2, nrow = 1)

# Make same histogram only for expressed genes
p1 <-ggplot(data = tfas_multicopy[tfas_multicopy$expressed_Tfas == "True",], aes(x = diff_longest_CDS)) +
  geom_histogram()
p2 <- ggplot(data = tlei_multicopy[tlei_multicopy$expressed_Tlei == "True",], aes(x = diff_longest_CDS)) +
  geom_histogram()
grid.arrange(p1, p2, nrow = 1)

# How many multicopy genes have both a start & stopcodon, either or none in Tfas
dim(tfas_multicopy[tfas_multicopy$startcodon == "True" & tfas_multicopy$stopcodon == "True",]) # 10,331 (86 % of multi copy genes)
dim(tfas_multicopy[(tfas_multicopy$startcodon == "True" & tfas_multicopy$stopcodon == "False") | 
                     (tfas_multicopy$startcodon == "False" & tfas_multicopy$stopcodon == "True"),]) # 1373 (11.4 % of multicopy genes)
dim(tfas_multicopy[tfas_multicopy$startcodon == "False" & tfas_multicopy$stopcodon == "False",]) # 310 genes (2.58 % of multicopy genes)
tfas_multicopy_incomplete <- tfas_multicopy[!(tfas_multicopy$startcodon == "True" & tfas_multicopy$stopcodon == "True"),]
write.table(tfas_multicopy_incomplete[,c(1,2)], file = "Incomplete_multicopy_genes_start-stopcodon_Tfas.txt", 
            quote = F, row.names = F)
# How many multicopy genes have both a start & stopcodon, either or none in Tlei
dim(tlei_multicopy[tlei_multicopy$startcodon == "True" & tlei_multicopy$stopcodon == "True",]) # 4011 (50 % of multicopy genes)
dim(tlei_multicopy[(tlei_multicopy$startcodon == "True" & tlei_multicopy$stopcodon == "False") | 
                     (tlei_multicopy$startcodon == "False" & tlei_multicopy$stopcodon == "True"),]) # 3382 (42 %)
dim(tlei_multicopy[tlei_multicopy$startcodon == "False" & tlei_multicopy$stopcodon == "False",]) # 618 (7.71 %)
tlei_multicopy_incomplete <- tlei_multicopy[!(tlei_multicopy$startcodon == "True" & tlei_multicopy$stopcodon == "True"),]
write.table(tlei_multicopy_incomplete[,c(1,2)], file = "Incomplete_multicopy_genes_start-stopcodon_Tlei.txt", 
            quote = F, row.names = F, sep = "\t")
# How many incomplete multicopy genes are expressed?
dim(tlei_multicopy_incomplete[tlei_multicopy_incomplete$expressed_Tfas == "True",]) # 2844 (71 % of incomplete genes)
dim(tlei_multicopy_incomplete[tlei_multicopy_incomplete$expressed_Tlei == "True",]) # 3014 (75 % of incomplete genes)

tfas_de_genes <- read.table("checklist_curated_orthologs_Tfas-Tlei.DEgenes_Tfas.txt", header = F)
colnames(tfas_de_genes) <- c("Gene_id",	"orthogroup",	"CDS_length",	"exon_count",	"startcodon",	"stopcodon",	
              "expressed_Tfas",	"expressed_Tlei",	"diff_largest_exon_count", "diff_longest_CDS")
dim(tfas_de_genes[tfas_de_genes$startcodon == "True" & tfas_de_genes$stopcodon == "True",]) # 725 (90.8 % of DE genes)
dim(tfas_de_genes[(tfas_de_genes$startcodon == "True" & tfas_de_genes$stopcodon == "False") | 
                     (tfas_de_genes$startcodon == "False" & tfas_de_genes$stopcodon == "True"),]) # 70 (8.77 % of DE genes)
dim(tfas_de_genes[tfas_de_genes$startcodon == "False" & tfas_de_genes$stopcodon == "False",]) # 3 genes (0.3 % of DE genes)
tfas_de_incomplete <- tfas_de_genes[!(tfas_de_genes$startcodon == "True" & tfas_de_genes$stopcodon == "True"),]
write.table(tfas_de_incomplete[,c(1,2)], file = "Incomplete_DE_genes_start-stopcodon_Tfas.txt", 
            quote = F, row.names = F, sep = "\t")

tlei_de_genes <- read.table("checklist_curated_orthologs_Tfas-Tlei.DEgenes_Tlei.txt", header = F)
colnames(tlei_de_genes) <- c("Gene_id",	"orthogroup",	"CDS_length",	"exon_count",	"startcodon",	"stopcodon",	
                             "expressed_Tfas",	"expressed_Tlei",	"diff_largest_exon_count", "diff_longest_CDS")
dim(tlei_de_genes[tlei_de_genes$startcodon == "True" & tlei_de_genes$stopcodon == "True",]) # 544 (75.3 % of DE genes)
dim(tlei_de_genes[(tlei_de_genes$startcodon == "True" & tlei_de_genes$stopcodon == "False") | 
                    (tlei_de_genes$startcodon == "False" & tlei_de_genes$stopcodon == "True"),]) # 169 (23.4 % of DE genes)
dim(tlei_de_genes[tlei_de_genes$startcodon == "False" & tlei_de_genes$stopcodon == "False",]) # 9 genes (1.2 % of DE genes)
tlei_de_incomplete <- tlei_de_genes[!(tlei_de_genes$startcodon == "True" & tlei_de_genes$stopcodon == "True"),]
write.table(tlei_de_incomplete[,c(1,2)], file = "Incomplete_DE_genes_start-stopcodon_Tlei.txt", 
            quote = F, row.names = F, sep = "\t")


### WITH NEW DATA FORMAT: looking at multiple stop codons and expression 
# on the exon level. Everything that is derived from "data" has expression
# values with relaxed filtering, any counts are seen as expression. Everything
# derived from "data2" has more stringent filtering and should be seen as the
# final result: here exons are only regarded as expressed when the mean CPM
# per sample is 0.001. CPM distributions were evaluated at the end of this file.
data <- read.delim("checklist_curated_orthologs_Tfas-Tlei.new062022.txt", header = T, sep = "\t")

data2 <- read.delim("checklist_curated_orthologs_Tfas-Tlei.new062022-2.txt",
                    header=T, sep="\t")
tfas_data <- data[startsWith(data$Gene_id, "Tfas"),]
tfas_data2 <- data2[startsWith(data2$Gene_id, "Tfas"),]
dim(tfas_data[tfas_data$startcodon == "True" & tfas_data$stopcodon == "one_stopcodon",]) # 23,618 (89.7 % of all genes are complete)
dim(tfas_data[tfas_data$stopcodon == "multiple_stopcodons",]) # Only one gene has multiple stopcodons. It is a unique gene of T. fasciculata
dim(tfas_data[tfas_data$stopcodon == "no_stopcodon",]) # 720 genes don't have a stopcodon
dim(tfas_data[tfas_data$startcodon == "False",]) # 2318 genes don't have a startcodon

tlei_data <- data[startsWith(data$Gene_id, "Tlei"),]
tlei_data2 <- data2[startsWith(data2$Gene_id, "Tlei"),]
dim(tlei_data[tlei_data$startcodon == "True" & tlei_data$stopcodon == "True",]) # 14,924 (63,3 % of all genes are complete)
dim(tlei_data[tlei_data$stopcodon == "multiple_stopcodons",]) # No genes with multiple stopcodons
dim(tlei_data[tlei_data$stopcodon == "no_stopcodon",]) # 6500 genes don't have a stopcodon
dim(tlei_data[tlei_data$startcodon == "False",]) # 2955 genes don't have a startcodon
dim(tlei_data[tlei_data$startcodon == "False" & tlei_data$stopcodon == "no_stopcodon",]) # 795 genes are fully incomplete

tlei_incomplete <- tlei_data[!(tlei_data$startcodon == "True" & tlei_data$stopcodon == "one_stopcodon"),]
write.table(tlei_incomplete, file = "Incomplete_genes_start-stopcodon_Tlei.txt", 
            quote = F, row.names = F, sep = "\t")
tlei_incomplete$expressed_Tlei <- as.numeric(tlei_incomplete$expressed_Tlei)
dim(tlei_incomplete[tlei_incomplete$expressed_Tlei == 1,]) # 8524 genes (98.4 % of all incomplete genes) have expression in all exons
tlei_incomplete2 <- tlei_data2[!(tlei_data2$startcodon == "True" & tlei_data2$stopcodon == "one_stopcodon"),]
dim(tlei_incomplete2[tlei_incomplete2$expressed_Tlei == "1.0",]) # After stronger filtering of expression, 5710 (65.9 %) of the incomplete genes are expressed across all exons.
dim(tlei_incomplete2[tlei_incomplete2$expressed_Tlei == "not_expressed",]) # After stronger filtering of expression, 2856 (32.9 %) of the incomplete genes are marked as not_expressed. This means that the majority of genes to filter out were not expressed or very lowly expressed across ALL exons.
dim(tlei_data2[tlei_data2$expressed_Tlei == "1.0",]) # 79.1 % of all Tlei genes are fully expressed in Tlei
dim(tlei_data2[tlei_data2$expressed_Tlei == "not_expressed",]) # 20 % of all Tlei genes are not at all expressed in Tlei

tfas_incomplete2 <- tfas_data2[!(tfas_data2$startcodon == "True" & tfas_data2$stopcodon == "one_stopcodon"),]
dim(tfas_incomplete2[tfas_incomplete2$expressed_Tfas == "1.0",]) # After stronger filtering of expression, 1798 (66.4 %) of the incomplete genes are expressed across all exons.
dim(tfas_incomplete2[tfas_incomplete2$expressed_Tfas == "not_expressed",]) # After stronger filtering of expression, 904 (33.4 %) of the incomplete genes are marked as not_expressed. This means that the majority of genes to filter out were not expressed or very lowly expressed across ALL exons.
dim(tfas_data2[tfas_data2$expressed_Tfas == "1.0",]) # 77.3 % of all tfas genes are fully expressed in tfas
dim(tfas_data2[tfas_data2$expressed_Tfas == "not_expressed",]) # 22.1 % of all Tlei genes are not at all expressed in Tlei




##############################
# Using EdgeR, I calculate CPM for the exons and find a threshold for
# further filtering of lowly expressed genes. This resulted in the checklist 
# 2, which then indicated which genes were incomplete and not expressed (less than 0.001 CPM)
counts_tlei <- read.delim("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/counts.Tfas_Tlei_6_timepoints.exons.ToTLEI.edited.txt", 
                          header = T, sep = "\t")
counts_tfas <- read.delim("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/mapped_to_Tfas/counts.Tfas_Tlei_6_timepoints.exons.edited.txt"
                          , header = T, sep = "\t")
library(edgeR)
library(stringr)
row.names(counts_tfas) <- counts_tfas[,1]
counts <- counts_tfas[,-1]
groups_list <- data.table::transpose(str_split(colnames(counts), "_"))[c(1,3,4)]
groups <- paste0(groups_list[[1]], "_", groups_list[[2]], "_", groups_list[[3]])
dyg<-DGEList(counts, group=groups)
dyg<-calcNormFactors(dyg, method="TMM")
normd <- cpm(dyg, normalized.lib.sizes = T)
normd_Tfas <- normd[,c(1:36)]
normd_Tlei <- normd[,c(37:72)]
# Here I explore the distribution of mean CPM per exon
x <- as.data.frame(rowMeans(normd_Tfas))
summary(x) # On average, an exon has 6 CPM, but the median is 0.6.
quantile(x$`rowMeans(normd_Tfas)`,probs = c(.1,.12,.5))
y <- as.data.frame(rowMeans(normd_Tlei))
summary(y) # On average, an exon has 6 CPM, but the median is 0.6.
quantile(y$`rowMeans(normd_Tlei)`,probs = c(.1,.14,.5))
# I only count as expressed the exons that have a mean CPM > 0.6. 
# I chose the median because that was a similar choice when filtering genes
write.table(normd, file = "/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/mapped_to_Tlei/counts.Tfas_Tlei_6_timepoints.exons.toTLEI.normalized-cpm.EdgeR.txt", 
            sep = "\t", quote = F, )

row.names(counts_tlei) <- counts_tlei[,1]
counts <- counts_tlei[,-1]
groups_list <- data.table::transpose(str_split(colnames(counts), "_"))[c(1,3,4)]
groups <- paste0(groups_list[[1]], "_", groups_list[[2]], "_", groups_list[[3]])
dyg<-DGEList(counts, group=groups)
dyg<-calcNormFactors(dyg, method="TMM")
normd <- cpm(dyg, normalized.lib.sizes = T)
normd_Tfas <- normd[,c(1:36)]
normd_Tlei <- normd[,c(37:72)]
# Here I explore the distribution of mean CPM per exon
x <- as.data.frame(rowMeans(normd_Tfas))
summary(x) # On average, an exon has 6 CPM, but the median is 0.6.
quantile(x$`rowMeans(normd_Tfas)`,probs = c(.1,.15,.5))
y <- as.data.frame(rowMeans(normd_Tlei))
summary(y) # On average, an exon has 6 CPM, but the median is 0.6.
quantile(y$`rowMeans(normd_Tlei)`,probs = c(.1,.15,.5))
# I only count as expressed the exons that have a mean CPM > 0.6. 
# I chose the median because that was a similar choice when filtering genes
write.table(normd, file = "/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/mapped_to_Tlei/counts.Tfas_Tlei_6_timepoints.exons.toTFAS.normalized-cpm.EdgeR.txt", 
            sep = "\t", quote = F, )


