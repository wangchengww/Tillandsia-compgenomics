#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2", "reshape2", "stringr","grid", "gridExtra")

# Load data
setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/")
counts <- read.table("counts.Tfas_Tlei_6_timepoints.normalized-cpm.EdgeR.logtransformed.txt", header = T, row.names = 1)
genes <- scan("Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF-cluster7.txt", character(), quote = "")
module <- str_split("Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF-cluster7.txt", "\\-|\\_|\\.")[[1]][11]

module_counts <- subset(counts, rownames(counts) %in% genes)
module_counts_Tfas <- module_counts[, c(1:36)]
module_counts_Tlei <- module_counts[, c(37:72)]
mod_size <- as.character(nrow(module_counts))

mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")

mean_count_Tfas$gene_id <- row.names(mean_count_Tfas)
mean_count_Tlei$gene_id <- row.names(mean_count_Tlei)

mean_count_Tfas_m <- melt(mean_count_Tfas, id.vars = "gene_id")
colnames(mean_count_Tfas_m) <- c("gene_id", "time", "count")
mean_count_Tlei_m <- melt(mean_count_Tlei, id.vars = "gene_id")
colnames(mean_count_Tlei_m) <- c("gene_id", "time", "count")

# Calculate total mean curve
means_Tfas <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
                         count = c(mean(mean_count_Tfas$`0100`), mean(mean_count_Tfas$`0500`),
                                   mean(mean_count_Tfas$`0900`), mean(mean_count_Tfas$`1300`),
                                   mean(mean_count_Tfas$`1700`), mean(mean_count_Tfas$`2100`)))

means_Tlei <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
                         count = c(mean(mean_count_Tlei$`0100`), mean(mean_count_Tlei$`0500`),
                                   mean(mean_count_Tlei$`0900`), mean(mean_count_Tlei$`1300`),
                                   mean(mean_count_Tlei$`1700`), mean(mean_count_Tlei$`2100`)))

# Highlight genes of interest
pdf(paste("Expression_curve_Tfas-vs-Tlei_TLEI-REF_", module, "_logtransformed.MEANS.pdf", sep = ""), width = 12, height = 8)
p1 <- ggplot(mean_count_Tfas_m, aes(x=time, y=count, group = gene_id)) +
  geom_line(data = means_Tfas, aes(group = 1), size = 1, color = "black") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(c(1.65,3.25))

p2 <- ggplot(mean_count_Tlei_m, aes(x=time, y=count, group = gene_id)) +
  geom_line(data = means_Tlei, aes(group = 1), size = 1, color = "black") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(c(4,5.6))

grid.arrange(p1, p2, nrow = 1, top = paste("Mean expression curve of ",module, " (", mod_size, ") in T. fasciculata
             (left) and T. leiboldiana (right)", sep = ""))
dev.off()
