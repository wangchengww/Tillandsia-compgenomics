#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2", "matrixStats", "stringr")

# Load arguments
# 1 is counts Tfas, 2 is counts Tlei, 3 is the gene of interest, 4 is orthology info.
args <- commandArgs(trailingOnly = TRUE)

# Load data
counts_Tfas <- read.table(args[[1]], header = T)
counts_Tlei <- read.table(args[[2]], header = T)
orthoinfo <- read.delim(args[[4]], sep = "\t", header = F)

gene = args[[3]]
func <- orthoinfo[orthoinfo$V1 == gene, 6]
func <- str_split(func, "=")[[1]][2]
orthology <- orthoinfo[orthoinfo$V1 == gene, c(8,9,10)]
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

pdf(paste("Expression_curve_", gene, ".pdf", sep = ""), width = 10, height = 8)
ggplot(mean_gene_expression, aes(x=time, y=mean_log_counts, group = species)) +
  geom_point(aes(color = species)) +
  geom_line(aes(color = species)) +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ggtitle(paste("Expression curve for ", gene, "\nFunction: " , func, "\nOrthology: ", orthology$V8, ":",
                orthology$V9, ":", orthology$V10, sep = "")) +
  geom_errorbar(aes(ymin=mean_log_counts-standard_dev, ymax=mean_log_counts+standard_dev), width=.2,
                position=position_dodge(0.05), colour = "darkgrey") +
  ylab("Read counts") +
  xlab ("Time")
dev.off()
