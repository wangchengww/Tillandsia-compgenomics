setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/duplication_stats/')
counts <- read.table("duplication_stats.per-orthogroup_EXONIC_absolute.txt", sep = ' ', header = T, row.names = 1)
counts <- as.data.frame(t(counts))
contingency <- as.data.frame(counts$all, row.names = rownames(counts))
contingency$non_dge <- counts$whole_genome - counts$all
colnames(contingency) <- c("dge", "non_dge")
chisq <- chisq.test(contingency)

library(corrplot)
corrplot(chisq$residuals, is.cor = FALSE)
contrib <- 100*chisq$residuals^2/chisq$statistic
round(contrib, 3)
