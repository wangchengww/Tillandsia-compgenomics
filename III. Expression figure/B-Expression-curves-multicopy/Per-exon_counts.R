setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/B-Expression-curves-multicopy/")

library(dplyr)
library(tidyselect)
library(ggplot2)

 dat <- read.table("counts.Tfas_Tlei_6_timepoints.exons.txt", head = T)

 pdf("Per-exon_counts_OG0000580_PFP.pdf")
gene1 <- subset(dat, startsWith(dat$Geneid, "Tfasc_v1.09221"))
gene1_tlei <- gene1[,43:78]
means <- rowMeans((gene1_tlei)) 
means <- as.data.frame(cbind(gene1$Geneid, as.numeric(means)))
means$pos <- c(1:19)
colnames(means) <- c("exon", "mean", "pos")
means$mean <- as.numeric(means$mean)
ggplot(means, aes(x=pos, y=mean)) +
  geom_bar(stat="identity") +
  ggtitle("Tfasc_v1.09221 (red) in Tlei")
summary(means)

gene1 <- subset(dat, startsWith(dat$Geneid, "Tfasc_v1.09221"))
gene1_tfas <- gene1[,7:42]
means <- rowMeans((gene1_tfas)) 
means <- as.data.frame(cbind(gene1$Geneid, as.numeric(means)))
means$pos <- c(1:19)
colnames(means) <- c("exon", "mean", "pos")
means$mean <- as.numeric(means$mean)
ggplot(means, aes(x=pos, y=mean)) +
  geom_bar(stat="identity") +
  ggtitle("Tfasc_v1.09221 (red) in Tfas")
summary(means)


gene1 <- subset(dat, startsWith(dat$Geneid, "Tfasc_v1.25891"))
gene1_tlei <- gene1[,43:78]
means <- rowMeans((gene1_tlei)) 
means <- as.data.frame(cbind(gene1$Geneid, as.numeric(means)))
means$pos <- c(1:19)
colnames(means) <- c("exon", "mean", "pos")
means$mean <- as.numeric(means$mean)
ggplot(means, aes(x=pos, y=mean)) +
  geom_bar(stat="identity") +
  ggtitle("Tfasc_v1.25891 (green) in Tlei")
summary(means)

gene1 <- subset(dat, startsWith(dat$Geneid, "Tfasc_v1.25891"))
gene1_tfas <- gene1[,7:42]
means <- rowMeans((gene1_tfas)) 
means <- as.data.frame(cbind(gene1$Geneid, as.numeric(means)))
means$pos <- c(1:19)
colnames(means) <- c("exon", "mean", "pos")
means$mean <- as.numeric(means$mean)
ggplot(means, aes(x=pos, y=mean)) +
  geom_bar(stat="identity") +
  ggtitle("Tfasc_v1.25891 (green) in Tfas")
summary(means)

gene1 <- subset(dat, startsWith(dat$Geneid, "Tfasc_v1.06440"))
gene1_tlei <- gene1[,43:78]
means <- rowMeans((gene1_tlei)) 
means <- as.data.frame(cbind(gene1$Geneid, as.numeric(means)))
means$pos <- c(1:19)
colnames(means) <- c("exon", "mean", "pos")
means$mean <- as.numeric(means$mean)
ggplot(means, aes(x=pos, y=mean)) +
  geom_bar(stat="identity") +
  ggtitle("Tfasc_v1.06440 (yellow) in Tlei")
summary(means)

gene1 <- subset(dat, startsWith(dat$Geneid, "Tfasc_v1.06440"))
gene1_tfas <- gene1[,7:42]
means <- rowMeans((gene1_tfas)) 
means <- as.data.frame(cbind(gene1$Geneid, as.numeric(means)))
means$pos <- c(1:19)
colnames(means) <- c("exon", "mean", "pos")
means$mean <- as.numeric(means$mean)
ggplot(means, aes(x=pos, y=mean)) +
  geom_bar(stat="identity") +
  ggtitle("Tfasc_v1.06440 (yellow) in Tfas")

dev.off()
summary(means)