setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/GO-term_Enrichment_Modules/')
setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/GO-term_Enrichment_Modules/')
library("ggplot2")
library(matrixStats)
library(stringr)
library(reshape2)

counts_Tfas <- read.table("../counts.Tfas.6_timepoints.filtr-normalized_DESEQ2.txt", header = T)
counts_Tlei <- read.table("../counts.Tlei.6_timepoints.filtr-normalized_DESEQ2.txt", header = T)
orthoinfo <- read.delim("../orthogroup_info_for_GOterm_enrichment.txt", sep = "\t", header = F)

gene = read.table("photorespiration_genes_unique_Tlei_signed18.txt", col.names = c("gene"))
gene_expression_Tfas <- subset(counts_Tfas, rownames(counts_Tfas) %in% gene$gene)
mean_gene_expression_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(gene_expression_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)]))                                           )
#sd_gene_expression_Tfas <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(gene_expression_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
#mean_gene_expression_Tfas <- cbind(mean_gene_expression_Tfas, sd_gene_expression_Tfas)
colnames(mean_gene_expression_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
mean_gene_expression_Tfas$Gene <- rownames(mean_gene_expression_Tfas)
mean_gene_expression_Tfas_m <- melt(mean_gene_expression_Tfas)
colnames(mean_gene_expression_Tfas_m) <- c("Gene","Time", "Counts")
mean_gene_expression_Tfas_m$Species <- "T.fasciculata"

gene_expression_Tlei <- subset(counts_Tlei, rownames(counts_Tlei) %in% gene$gene)
mean_gene_expression_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(gene_expression_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
#sd_gene_expression_Tlei <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(gene_expression_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
#mean_gene_expression_Tlei <- cbind(mean_gene_expression_Tlei, sd_gene_expression_Tlei)
colnames(mean_gene_expression_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")
mean_gene_expression_Tlei$Gene <- rownames(mean_gene_expression_Tlei)
mean_gene_expression_Tlei_m <- melt(mean_gene_expression_Tlei)
colnames(mean_gene_expression_Tlei_m) <- c("Gene", "Time", "Counts")
mean_gene_expression_Tlei_m$Species <- "T.leiboldiana"
mean_gene_expression <- rbind(mean_gene_expression_Tfas_m, mean_gene_expression_Tlei_m)
mean_gene_expression$Counts <- log(mean_gene_expression$Counts+1)

ggplot(mean_gene_expression, aes(x=Time, y=Counts, group = interaction(Gene, Species))) + 
  geom_point(aes(color = Species, shape=Gene)) +
  geom_line(aes(color = Species)) +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") 
  ggtitle() 
  geom_errorbar(aes(ymin=mean_log_counts-standard_dev, ymax=mean_log_counts+standard_dev), width=.2,
                position=position_dodge(0.05), colour = "darkgrey")



