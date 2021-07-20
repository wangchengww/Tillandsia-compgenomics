# Studying differences in gene family sizes between Tfas and Tlei

library(ggplot2)

setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
counts <- read.table("orthogroups_Tfas_Tlei_Acom.counts.no_TEs.size_corrections.txt", sep = '\t')
colnames(counts) <- c("og_id", "Acom", "Tfas", "Tlei")

# Filter out unique orthogroups
counts_Tfas_Tlei <- counts[counts$Tfas != 0 & counts$Tlei != 0,]
counts_Tfas_Tlei$diff <- counts_Tfas_Tlei$Tfas - counts_Tfas_Tlei$Tlei


ggplot(counts_Tfas_Tlei, aes(x=Tfas, y=Tlei)) + geom_point(size = .5)
