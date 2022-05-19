library(tidyverse)
library(rstatix)
library(ggpubr)
library(reshape2)

setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/")
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/")
counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.logtransformed.txt", header = T, row.names = 1)
genes <- scan("Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic-cluster7.txt", character(), quote = "")
module <- str_split("Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic-cluster7.txt", "\\-|\\_|\\.")[[1]][12]

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

# Calculate average differences per group
mean_count_Tfas_m$species <- "Tfas"
mean_count_Tlei_m$species <- "Tlei"
test <- rbind(mean_count_Tfas_m, mean_count_Tlei_m)
# This shows that none of the counts are normally distributed - we will have to perform a 
# Mann-Whitney U test
for (i in c(1:6)){
  t = as.vector(unique(test$time))[i]
  x <- test[test$time == t,]
  s <- x %>%
    group_by(species) %>%
    shapiro_test(count)
  ggqqplot(x, x = "count", facet.by = "species")
  name <- paste0("shapiro_",t)
  assign(name, s)
}

cluster_stats <- data.frame()
for (i in c(1:6)){
  t = as.vector(unique(test$time))[i]
  x <- test[test$time == t,]
  mean_Tfas <- round(mean(x[x$species == "Tfas",]$count), 2)
  mean_Tlei <- round(mean(x[x$species == "Tlei",]$count), 2)
  wilcox <- wilcox.test(count ~ species, data=x)
  line <- c(t, mean_Tfas, mean_Tlei, wilcox$p.value)
  cluster_stats <- rbind(cluster_stats, line)
  colnames(cluster_stats) <- c("Time", "Mean_log(CPM)_Tfas", "Mean_log(CPM)_Tlei", "p-value")
}
cluster_stats$Difference <- as.numeric(cluster_stats$`Mean_log(CPM)_Tfas`) - as.numeric(cluster_stats$`Mean_log(CPM)_Tlei`)
mean(cluster_stats$Difference)
write.table(cluster_stats, file = paste0("Mean_Expression_per_Timepoint_", module, ".EXONIC.txt"),
            quote = F, sep = "\t")
            