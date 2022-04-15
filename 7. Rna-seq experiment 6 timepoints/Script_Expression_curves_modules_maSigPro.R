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
counts <- read.table(args[[1]], header = T, row.names = 1)
genes <- scan(args[[2]], character(), quote = "")
module <- str_split(args[[2]], "\\-|\\_|\\.")[[1]][12]
goterms <- read.delim(args[[3]], sep = "\t", header = T)

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
write.table(mean_count_Tfas_m, file = paste0("Mean_Expression_counts_T.fasciculata-", module, ".txt"))
write.table(mean_count_Tlei_m, file = paste0("Mean_Expression_counts_T.leiboldiana-", module, ".txt"))

malate <- goterms[grep("malate", goterms$Term), 6]
malate <- unlist(str_split(malate, ", "))
malate_expr_Tfas <- subset(mean_count_Tfas_m, gene_id %in% malate)
malate_expr_Tlei <- subset(mean_count_Tlei_m, gene_id %in% malate)

phospho <- goterms[grep("phosphoenolpyruvate", goterms$Term), 6]
phospho <- unlist(str_split(phospho, ", "))
phospho_expr_Tfas <- subset(mean_count_Tfas_m, gene_id %in% phospho)
phospho_expr_Tlei <- subset(mean_count_Tlei_m, gene_id %in% phospho)

vacuole <- goterms[grep("vacuolar acidification", goterms$Term), 6]
vacuole <- unlist(str_split(vacuole, ", "))
vacuole_expr_Tfas <- subset(mean_count_Tfas_m, gene_id %in% vacuole)
vacuole_expr_Tlei <- subset(mean_count_Tlei_m, gene_id %in% vacuole)

stomata <- goterms[grep("circadian", goterms$Term), 6]
stomata <- unlist(str_split(stomata, ", "))
stomata_expr_Tfas <- subset(mean_count_Tfas_m, gene_id %in% stomata)
stomata_expr_Tlei <- subset(mean_count_Tlei_m, gene_id %in% stomata)

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
pdf(paste("Expression_curve_Tfas-vs-Tlei_TLEI-REF-exonic_", module, "_logtransformed.pdf", sep = ""), width = 12, height = 8)
p1 <- ggplot(mean_count_Tfas_m, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means_Tfas, aes(group = 1), size = 1, color = "black") +
  geom_line(data = malate_expr_Tfas, aes(group = gene_id), size = 1, color = "darkgreen") +
  geom_line(data = phospho_expr_Tfas, aes(group = gene_id), size = 1, color = "blue") +
  geom_line(data = vacuole_expr_Tfas, aes(group = gene_id), size = 1, color = "purple") +
  geom_line(data = stomata_expr_Tfas, aes(group = gene_id), size = 1, color = "red") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(c(ymin, ymax))

p2 <- ggplot(mean_count_Tlei_m, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means_Tlei, aes(group = 1), size = 1, color = "black") +
  geom_line(data = malate_expr_Tlei, aes(group = gene_id), size = 1, color = "darkgreen") +
  geom_line(data = phospho_expr_Tlei, aes(group = gene_id), size = 1, color = "blue") +
  geom_line(data = vacuole_expr_Tlei, aes(group = gene_id), size = 1, color = "purple") +
  geom_line(data = stomata_expr_Tlei, aes(group = gene_id), size = 1, color = "red") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(ymin, ymax)

grid.arrange(p1, p2, nrow = 1, top = paste("Expression curve of ",module, " (", mod_size, ") in T. fasciculata
             (left) and T. leiboldiana (right)", sep = ""))
dev.off()
