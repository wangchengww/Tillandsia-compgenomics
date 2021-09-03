# Studying differences in gene family sizes between Tfas and Tlei

library(ggplot2)

setwd("Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
counts <- read.table("orthogroups_Tfas_Tlei_Acom.counts.no_TEs.size_corrections.no_plastid-mito-ribo2.txt", sep = '\t')
colnames(counts) <- c("og_id", "Acom", "Tfas", "Tlei")

# Filter out unique orthogroups
counts_Tfas_Tlei <- counts[counts$Tfas != 0 & counts$Tlei != 0,]
ggplot(counts_Tfas_Tlei, aes(x=Tfas, y=Tlei)) + geom_point(size = .5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Per-species gene counts in multi-copy orthogroups") +
  ylab(label = "T. leiboldiana") +
  xlab(label = "T. fasciculata") +
  theme_bw()

# Add colour gradient to show number of datapoints with the same count combination
counts_Tfas_Tlei_multi <- counts_Tfas_Tlei[!(counts_Tfas_Tlei$Tfas == 1 & counts_Tfas_Tlei$Tlei == 1),]
ggplot(counts_Tfas_Tlei_multi) + geom_hex(aes(Tfas, Tlei), bins = 100) +
  labs(title = "Per-species gene counts in multi-copy orthogroups") +
  ylab(label = "T. leiboldiana") +
  xlab(label = "T. fasciculata") +
  scale_fill_continuous(type = "viridis") +
  theme_bw()
  theme_bw()

write.table(counts_Tfas_Tlei_npl_multi, file = "orthogroup_selection_for_GO_term_all.txt", sep = "\t")


