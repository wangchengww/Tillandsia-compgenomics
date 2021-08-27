# Studying differences in gene family sizes between Tfas and Tlei

library(ggplot2)

setwd("Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
counts <- read.table("orthogroups_Tfas_Tlei_Acom.counts.no_TEs.size_corrections.no_plastid-mito-ribo.txt", sep = '\t')
colnames(counts) <- c("og_id", "Acom", "Tfas", "Tlei")

# Filter out unique orthogroups
counts_Tfas_Tlei <- counts[counts$Tfas != 0 & counts$Tlei != 0,]
ggplot(counts_Tfas_Tlei, aes(x=Tfas, y=Tlei)) + geom_point(size = .5)

# Same without plastid annotations
counts_noplastid <- read.table("orthogroups_Tfas_Tlei_Acom.counts.no_TEs.size_corrections-bothspecies.no_plastid-mito-ribo.txt",
                               sep = "\t")
colnames(counts_noplastid) <- c("og_id", "Acom", "Tfas", "Tlei")
counts_Tfas_Tlei_npl <- counts_noplastid[counts_noplastid$Tfas != 0 & counts_noplastid$Tlei != 0,]
counts_Tfas_Tlei_npl <- counts_Tfas_Tlei_npl[counts_Tfas_Tlei_npl$og_id != "OG0000015",]

ggplot(counts_Tfas_Tlei_npl, aes(x=Tfas, y=Tlei)) + geom_point(size = .5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Per-species gene counts in multi-copy orthogroups") +
  ylab(label = "T. leiboldiana") +
  xlab(label = "T. fasciculata") +
  theme_bw()

counts_Tfas_Tlei_npl_multi <- counts_Tfas_Tlei_npl[!(counts_Tfas_Tlei_npl$Tfas == 1 & counts_Tfas_Tlei_npl$Tlei == 1),]
write.table(counts_Tfas_Tlei_npl_multi, file = "orthogroup_selection_for_GO_term_all.txt", sep = "\t")

chisq.test(counts_Tfas_Tlei_npl$Tfas, counts_Tfas_Tlei_npl$Tlei, correct=FALSE)
  
