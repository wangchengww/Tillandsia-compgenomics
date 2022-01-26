setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
library("ggplot2")
library("reshape2")

counts <- read.table("counts.Tfas.6_timepoints.filtr-normalized_DESEQ2.txt", header = T, row.names = 1)
genes <- scan("Genes-for-Enrichment_T.fasciculata_signed18_vst10c4s_lightcyan.txt", character(), quote = "")

module_counts<- subset(counts, rownames(counts) %in% genes)

mean_count <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count) <- c("0100", "0500", "0900","1300", "1700", "2100")
mean_count$gene_id <- row.names(mean_count)
mean_count_m <- melt(mean_count, id.vars = "gene_id")
colnames(mean_count_m) <- c("gene_id", "time", "count")

ggplot(mean_count_m, aes(x=time, y=count)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.3)
