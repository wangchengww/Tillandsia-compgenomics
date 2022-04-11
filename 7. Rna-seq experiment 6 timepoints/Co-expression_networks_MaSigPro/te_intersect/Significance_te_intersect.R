library(tidyverse)
library(rstatix)

setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/te_intersect/')

counts <- read.table('Tillandsia_fasciculata_GENE-TE-intersection.counts.txt', sep = '\t')
colnames(counts) <- c('chr', 'start', 'end', 'dir','name', 'te_counts')
genes <- scan('Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic.txt', character())
DE_counts <- subset(counts, counts$name %in% genes)
counts$type <- ifelse(counts$name %in% genes, "subset", "total")

# Calculate whether differentially expressed genes have on average higher rates of intronic TEs than
# the background using the welch test
stat.test <- counts %>% 
  t_test(te_counts ~ type) %>%
  add_significance()
# p-value = 0.00103

# Calculate whether the percentage of genes with a TE insertion is significantly higher in DE genes
contingency <- data.frame(rbind(c(sum(counts$te_counts == 0), 
                                    sum(counts[counts$type == "subset",]$te_counts == 0)),
                                  c(sum(counts$te_counts != 0),
                                    sum(counts[counts$type == "subset",]$te_counts != 0))))
colnames(contingency) <- c("total", "subset")
rownames(contingency) <- c("no_te", "te")

chisq <- chisq.test(contingency)

