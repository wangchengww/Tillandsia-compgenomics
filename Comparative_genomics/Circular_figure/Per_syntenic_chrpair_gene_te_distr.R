### Comparison of syntenic chromosomes in Tfas and Tlei and the positions of gene-rich and TE-rich 
#### areas

## The idea here is to create per-chromosome pair flat plots to see more clearly whether gene-rich 
## areas have shifted between both species

setwd('/Users/clara/bio-info_phd/Comparative_genomics/Circular_figure/')

# Read in genic content
Tfas_onetoone <- read.table("Tfas_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
Tlei_onetoone <- read.table("Tlei_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
Tfas_onetoone$V2 <- (Tfas_onetoone$V2)/1000000
Tlei_onetoone$start_window <- (Tlei_onetoone$start_window)/1000000

# Read in repetitive content
Tfas_rep <- read.table("TE_content_Tfas_per1MBwindow_py.txt", header = T)
Tfas_rep$start_window <- (Tfas_rep$start_window)/1000000
Tlei_rep <- read.table("TE_content_Tlei_per1MBwindow_py.txt", header = T)
Tlei_rep$start_window <- (Tlei_rep$start_window)/1000000
summary(Tfas_rep)
# Chr1 Tfas
# Genic content
chr1 <- Tfas_onetoone[Tfas_onetoone$chrom == "Scaffold_2199",]
plot(chr1$V4~chr1$V2, type="l", lwd=2.5,lty=1, ylim=c(0,100), 
     xlab="position (MB)",
     ylab="Gene count", 
     xaxs = "i",
     yaxs = "i",
     col="seagreen",xaxt="n")+
  title(main = "Gene count per 1 MB window on Tfas_chr1",)
axis(1, xaxp  = c(0, 32, 8), las=1)
axis(1, tck=-0.015, col.ticks="black", xaxp = c(0,32, 32), labels = F)
polygon(c(chr1$V2[chr1$V2>=0], max(chr1$V2), 0), c(chr1$V4[chr1$V2>=0], 0, 0), 
        col=adjustcolor("seagreen",alpha.f=0.7))
grid(nx = 8, ny = 5)

# Repetitive content
Tfas_chr1_rep <- Tfas_rep[Tfas_rep$chrom == "Tfas_chr1",]
plot(Tfas_chr1_rep$perc_te~Tfas_chr1_rep$start_window, type="l", lwd=2.5,lty=1, ylim=c(0,100), 
     xlab="position (MB)",
     ylab="Proportion of masked bases", 
     xaxs = "i",
     yaxs = "i",
     col="orange",xaxt="n")+
  title(main = "Repetitive content per 1 MB window on Tfas_chr1",)
axis(1, xaxp  = c(0, 32, 8), las=1)
axis(1, tck=-0.015, col.ticks="black", xaxp = c(0,32, 32), labels = F)
polygon(c(Tfas_chr1_rep$start_window[Tfas_chr1_rep$start_window>=0], max(Tfas_chr1_rep$start_window), 0), c(Tfas_chr1_rep$perc_te[Tfas_chr1_rep$start_window>=0], 0, 0), 
        col=adjustcolor("orange",alpha.f=0.7))
grid(nx = 8, ny = 5)


# Chr2 Tlei
# Genic content
chr2_lei <- Tlei_onetoone[Tlei_onetoone$chrom == "Scaffold_10370",]
plot(chr2_lei$gene_count~chr2_lei$start_window, type="l", lwd=2.5,lty=1, ylim=c(0,100), 
     xlab="position (MB)",
     ylab="Gene count", 
     xaxs = "i",
     yaxs = "i",
     col="seagreen",xaxt="n")+
  title(main = "Gene count per 1 MB window on Tlei_chr2",)
axis(1, xaxp  = c(0, 60, 12), las=1)
axis(1, tck=-0.015, col.ticks="black", xaxp = c(0,58, 58), labels = F)
polygon(c(chr2_lei$start_window[chr2_lei$start_window>=0], max(chr2_lei$start_window), 0), c(chr2_lei$gene_count[chr2_lei$start_window>=0], 0, 0), 
        col=adjustcolor("seagreen",alpha.f=0.7))
grid()

# Repetitive content
Tlei_chr2_rep <- Tlei_rep[Tlei_rep$chrom == "Tlei_chr2",]
plot(Tlei_chr2_rep$perc_te~Tlei_chr2_rep$start_window, type="l", lwd=2.5,lty=1, ylim=c(0,100), 
     xlab="position (MB)",
     ylab="Proportion of masked bases", 
     xaxs = "i",
     yaxs = "i",
     col="orange",xaxt="n")+
  title(main = "Repetitive content per 1 MB window on Tlei_chr2",)
axis(1, xaxp  = c(0, 60, 12), las=1)
axis(1, tck=-0.015, col.ticks="black", xaxp = c(0,58, 58), labels = F)
polygon(c(Tlei_chr2_rep$start_window[Tlei_chr2_rep$start_window>=0], max(Tlei_chr2_rep$start_window), 0), c(Tlei_chr2_rep$perc[Tlei_chr2_rep$start_window>=0], 0, 0), 
        col=adjustcolor("orange",alpha.f=0.7))
grid()

dev.off()

#---
#Chr26 Tlei

chr21_lei <- Tlei_onetoone[Tlei_onetoone$chrom == "Scaffold_10430",]
plot(chr21_lei$gene_count~chr21_lei$start_window, type="l", lwd=2.5,lty=1, ylim=c(0,100), 
     xlab="position (MB)",
     ylab="Gene count", 
     xaxs = "i",
     yaxs = "i",
     col="seagreen",xaxt="n")+
  title(main = "Gene count per 1 MB window on Tlei_chr21",)
axis(1, xaxp  = c(0, 12, 6), las=1)
axis(1, tck=-0.015, col.ticks="black", xaxp = c(0,58, 58), labels = F)
axis(1, tck=-0.015, col.ticks="black", xaxp = c(0,58, 58), labels = F)
polygon(c(chr21_lei$start_window[chr21_lei$start_window>=0], max(chr21_lei$start_window), 0), c(chr21_lei$gene_count[chr21_lei$start_window>=0], 0, 0), 
        col=adjustcolor("seagreen",alpha.f=0.7))
grid()

# Repetitive content
Tlei_chr21_rep <- Tlei_rep[Tlei_rep$chrom == "Tlei_chr21",]
plot(Tlei_chr21_rep$perc_te~Tlei_chr21_rep$start_window, type="l", lwd=2.5,lty=1, ylim=c(0,100), 
     xlab="position (MB)",
     ylab="Proportion of masked bases", 
     xaxs = "i",
     yaxs = "i",
     col="orange",xaxt="n")+
  title(main = "Repetitive content per 1 MB window on Tlei_chr21",)
axis(1, xaxp  = c(0, 10, 5), las=1)
axis(1, tck=-0.015, col.ticks="black", xaxp = c(0,58, 58), labels = F)
polygon(c(Tlei_chr21_rep$start_window[Tlei_chr21_rep$start_window>=0], max(Tlei_chr21_rep$start_window), 0), c(Tlei_chr21_rep$perc[Tlei_chr21_rep$start_window>=0], 0, 0), 
        col=adjustcolor("orange",alpha.f=0.7))
grid()