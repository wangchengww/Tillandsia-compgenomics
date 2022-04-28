#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2", "reshape2", "stringr","grid", "gridExtra")

# Load arguments
# 1 is counts, 2 is the subset of genes of each module, 3 is the GOterms
args <- commandArgs(trailingOnly = TRUE)

# Load data
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/B-Expression-curves-multicopy/')
counts <- read.table(args[[1]], header = T, row.names = 1)
counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.logtransformed.txt", header = T, row.names = 1)
genes <- scan(args[[2]], character(), quote = "")
genes <- scan("OG0001504_PEPC_gene_list.txt", character(), quote = "")
module <- str_split("OG0001504_PEPC_gene_list.txt", "_")[[1]][1]

module_counts <- subset(counts, rownames(counts) %in% genes)
module_counts_Tfas <- module_counts[, c(1:36)]
module_counts_Tlei <- module_counts[, c(37:72)]
mod_size <- as.character(nrow(module_counts))

# Melt the count data
module_counts_Tfas$gene_id <- row.names(module_counts_Tfas)
module_counts_Tfas_m <- melt(module_counts_Tfas, id.vars = "gene_id")
module_counts_Tfas_m$time <- c(rep("0100", mod_size), rep("0500",mod_size), rep("0900",mod_size),
                               rep("1300",mod_size), rep("1700",mod_size), rep("2100",mod_size))
module_counts_Tfas_m$sample <- c(rep("A", as.numeric(mod_size)*6), rep("B",as.numeric(mod_size)*6),
                                 rep("C",as.numeric(mod_size)*6), rep("D",as.numeric(mod_size)*6),
                                 rep("E",as.numeric(mod_size)*6), rep("F",as.numeric(mod_size)*6))
colnames(module_counts_Tfas_m) <- c("gene_id", "id", "count", "time", "sample")

module_counts_Tlei$gene_id <- row.names(module_counts_Tlei)
module_counts_Tlei_m <- melt(module_counts_Tlei, id.vars = "gene_id")
module_counts_Tlei_m$time <- c(rep("0100", mod_size), rep("0500",mod_size), rep("0900",mod_size),
                               rep("1300",mod_size), rep("1700",mod_size), rep("2100",mod_size))
module_counts_Tlei_m$sample <- c(rep("A", as.numeric(mod_size)*6), rep("B",as.numeric(mod_size)*6),
                                 rep("C",as.numeric(mod_size)*6), rep("D",as.numeric(mod_size)*6),
                                 rep("E",as.numeric(mod_size)*6), rep("F",as.numeric(mod_size)*6))
colnames(module_counts_Tlei_m) <- c("gene_id", "id", "count", "time", "sample")

# Calculate mean counts per time point
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

# Calculate scale
max_Tfas <- max(mean_count_Tfas_m$count)
min_Tfas <- min(mean_count_Tfas_m$count)
max_Tlei <- max(mean_count_Tlei_m$count)
min_Tlei <- min(mean_count_Tlei_m$count)

if (max_Tfas > max_Tlei){
  ymax = max_Tfas + 1
} else {
  ymax = max_Tlei + 1
}

if (min_Tfas < min_Tlei){
  ymin = min_Tfas - 1
} else {
  ymin = min_Tlei - 1
}

# Highlight genes of interest
pdf(paste("Expression_curve_", module, "_exonic.pdf", sep = ""), width = 12, height = 8)
p1 <- ggplot(mean_count_Tfas_m, aes(x=time, y=count, group = gene_id)) +
  geom_point(aes(color=gene_id)) +
  geom_line(aes(color=gene_id)) +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  scale_color_manual(values=c('purple','purple4')) +
  ylim(c(ymin, ymax)) +
  ylab("Median-centered log(CPM)") +
  xlab("Time") +
  labs(fill='Gene ID')

p2 <- ggplot(mean_count_Tlei_m, aes(x=time, y=count, group = gene_id)) +
  geom_point(aes(color=gene_id)) +
  geom_line(aes(color=gene_id)) +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  scale_color_manual(values=c('purple','purple4')) +
  ylim(ymin, ymax) +
  ylab("Median-centered log(CPM)") +
  xlab("Time") +
  labs(fill='Gene ID')

grid.arrange(p1, p2, nrow = 1, top = paste("Expression curve of ",module, " (", mod_size, ") in T. fasciculata
             (left) and T. leiboldiana (right)", sep = ""))
dev.off()
