setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/GO-term_Enrichment_Modules/')
library("ggplot2")
library(matrixStats)

counts_Tfas <- read.table("../counts.Tfas.6_timepoints.filtr-normalized_DESEQ2.txt", header = T)
counts_Tlei <- read.table("../counts.Tlei.6_timepoints.filtr-normalized_DESEQ2.txt", header = T)
gene = "Tfasc_v1.09872"
gene_expression_Tfas <- counts_Tfas[rownames(counts_Tfas) == gene,]
mean_gene_expression_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(gene_expression_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)]))
                                           )
sd_gene_expression_Tfas <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(gene_expression_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
mean_gene_expression_Tfas <- cbind(mean_gene_expression_Tfas, sd_gene_expression_Tfas)
mean_gene_expression_Tfas$time <- c("0100", "0500", "0900","1300", "1700", "2100")
colnames(mean_gene_expression_Tfas) <- c("mean_log_counts", "standard_dev", "time")
mean_gene_expression_Tfas$species <- "T.fasciculata"

gene_expression_Tlei <- counts_Tlei[rownames(counts_Tlei) == gene,]
mean_gene_expression_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(gene_expression_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
sd_gene_expression_Tlei <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(gene_expression_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
mean_gene_expression_Tlei <- cbind(mean_gene_expression_Tlei, sd_gene_expression_Tlei)
mean_gene_expression_Tlei$time <- c("0100", "0500", "0900","1300", "1700", "2100")
colnames(mean_gene_expression_Tlei) <- c("mean_log_counts", "standard_dev", "time")
mean_gene_expression_Tlei$species <- "T.leiboldiana"
mean_gene_expression <- rbind(mean_gene_expression_Tfas, mean_gene_expression_Tlei)

ggplot(mean_gene_expression, aes(x=time, y=mean_log_counts, group = species)) + 
  geom_point(aes(color = species)) +
  geom_line(aes(color = species)) +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ggtitle(paste("Expression curve for", gene, "- function: PEP Carboxylase")) +
  geom_errorbar(aes(ymin=mean_log_counts-standard_dev, ymax=mean_log_counts+standard_dev), width=.2,
                position=position_dodge(0.05))
  ylim(c(0,270000))


