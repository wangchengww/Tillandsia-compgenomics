### Comparison of syntenic chromosomes in Tfas and Tlei and the positions of gene-rich and TE-rich 
#### areas

## The idea here is to create per-chromosome pair flat plots to see more clearly whether gene-rich 
## areas have shifted between both species

setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure/')
library("grid")
library("gridExtra")
# Read in genic content
gene_content <- read.table("Gene_counts_per_1MB_windows.Tfas-Tlei.mainScaffolds.curatedOGs.txt", header = T)
gene_content$start_window <- (gene_content$start_window)/1000000

# Read in repetitive content
rep_content <- read.table("TE_content_Tfas-Tlei_per1MB-window_python.txt", header = T)
rep_content$start_window <- (rep_content$start_window)/1000000

# Chr18 Tfas
# Genic content
chr18 <- gene_content[gene_content$chrom == "Tfas_chr18",]
par(mfrow=c(1,1))
plot(chr18$gene_counts~chr18$start_window, type="l", lwd=2.5,lty=1, ylim=c(0,130), 
     xlab="position (MB)",
     ylab="Gene count", 
     xaxs = "i",
     yaxs = "i",
     col="seagreen",xaxt="n")+
  title(main = "Gene count per 1 MB window on Tfas_chr18",)  
axis(1, xaxp  = c(0, 32, 8), las=1)
axis(1, tck=-0.015, col.ticks="black", xaxp = c(0,32, 32), labels = F)
polygon(c(chr18$start_window[chr18$start_window>=0], max(chr18$start_window), 0), c(chr18$gene_counts[chr18$start_window>=0], 0, 0), 
        col=adjustcolor("seagreen",alpha.f=0.7))
grid(nx = 8, ny = 5)
rect(xleft= 5.6,xright = 11,,ybottom=0,ytop=130, density=10, col = "blue")
rect(xleft= 17,xright = 21,,ybottom=0,ytop=130, density=10, col = "blue")
# Repetitive content
Tfas_chr18_rep <- rep_content[rep_content$chrom == "Tfas_chr18",]
 plot(Tfas_chr18_rep$perc_te~Tfas_chr18_rep$start_window, type="l", lwd=2.5,lty=1, ylim=c(0,100), 
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
rect(xleft= 5.6,xright = 11,,ybottom=0,ytop=130, density=10, col = "blue")
rect(xleft= 17,xright = 21,,ybottom=0,ytop=130, density=10, col = "blue")

grid.arrange(p1, p2, nrow = 1)
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