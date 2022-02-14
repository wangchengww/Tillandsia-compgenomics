#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2", "reshape2", "stringr")

# Load arguments
# 1 is counts, 2 is the subset of genes of each module, 3 is the GOterms
args <- commandArgs(trailingOnly = TRUE)

# Load data
counts <- read.table(args[[1]], header = T, row.names = 1)
genes <- scan(args[[2]], character(), quote = "")
module <- str_split(args[[2]], "\\_|\\.")[[1]][5]
species <- str_split(args[[2]], "\\_|\\.")[[1]][3]
network <- str_split(args[[2]], "\\_|\\.")[[1]][4]
goterms <- read.delim(args[[3]], sep = "\t", header = T)

# Subset counts
module_counts <- subset(counts, rownames(counts) %in% genes)
mod_size <- as.character(nrow(module_counts))

# Calculate mean of each gene across individuals for each timepoint
mean_count <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count) <- c("0100", "0500", "0900","1300", "1700", "2100")
mean_count$gene_id <- row.names(mean_count)
mean_count_m <- melt(mean_count, id.vars = "gene_id")
colnames(mean_count_m) <- c("gene_id", "time", "count")


malate <- goterms[grep("malate", goterms$Term), 6]
malate <- unlist(str_split(malate, ", "))
malate_expr <- subset(mean_count_m, gene_id %in% malate)

phospho <- goterms[grep("phosphoenolpyruvate", goterms$Term), 6]
phospho <- unlist(str_split(phospho, ", "))
phospho_expr <- subset(mean_count_m, gene_id %in% phospho)

vacuole <- goterms[grep("vacuolar acidification", goterms$Term), 6]
vacuole <- unlist(str_split(vacuole, ", "))
vacuole_expr <- subset(mean_count_m, gene_id %in% vacuole)

stomata <- goterms[grep("stomata", goterms$Term), 6]
stomata <- unlist(str_split(stomata, ", "))
stomata_expr <- subset(mean_count_m, gene_id %in% stomata)

# Calculate total mean curve
means <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
               count = c(mean(mean_count$`0100`), mean(mean_count$`0500`),
                 mean(mean_count$`0900`), mean(mean_count$`1300`),
                 mean(mean_count$`1700`), mean(mean_count$`2100`)))


# Highlight genes of interest
pdf(paste("Expression_curve_", species, "_", network, "_", module, "_logtransformed.pdf", sep = ""), width = 8, height = 10)
cols <- c('Malate'='red', 'PhosphoEnolPyruvate'='blue', 'Vacuole'='purple', 'Stomata'='darkgreen')
ggplot(mean_count_m, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means, aes(group = 1), size = 1, color = "black") +
  geom_line(data = malate_expr, aes(group = gene_id), size = 1, color = "red") +
  geom_line(data = phospho_expr, aes(group = gene_id), size = 1, color = "blue") +
  geom_line(data = vacuole_expr, aes(group = gene_id), size = 1, color = "purple") +
  geom_line(data = stomata_expr, aes(group = gene_id), size = 1, color = "darkgreen") +
  ggtitle(paste("Expression curve of ",module, " (", mod_size, ") in ", species, ", ", network, sep = "")) +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed")
dev.off()
