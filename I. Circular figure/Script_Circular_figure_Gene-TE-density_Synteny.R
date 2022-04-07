#!/usr/bin/Rscript --vanilla

# Code for large circular figure - Genome Paper T.fasciculata / T. leiboldiana
# This Rcode creates a circular figure with the 25 / 26 main scaffolds of both assemblies as halves of the circle. The circle has 3 tracks, a Gene density, a TE density and a synteny plot

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("circlize", "stringr", "RColorBrewer")

setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure")

# Load arguments
# 1 is counts, 2 is the subset of genes of each module, 3 is the GOterms
args <- commandArgs(trailingOnly = TRUE)
output_name <- args[[5]]
# Read in the complete list of chromosomes with all needed info:
#chrom <- read.table("chromosomes_Tfas_Tlei.coordinates-circle.txt", header = T, sep = "\t")
chrom <- read.table(args[[1]], header = T, sep = "\t")

# Make matrix of start and end position to initialize the circular plot
pos <- cbind(rep(0, 51), chrom$size)

# Make plot
pdf(paste0(output_name, ".pdf"))
print("Initializing the plot...")
circos.clear()
# Plot initialization
circos.par("track.height"=0.05, cell.padding=c(0.02, 0, 0.02, 0),
           start.degree=90,gap.degree=1,points.overflow.warning=FALSE,
           track.margin=c(0,0))
circos.initialize(sectors = chrom$name, xlim = pos)

# Add blocks representing the chromosomes
circos.track(chrom$name, ylim = c(0, 1))

# Add chromosome names
is.even <- function(x) x %% 2 == 0
n = 0
m = 0
for (i in 1:nrow(chrom)){
  name=chrom[i,1]
  if (is.even(i) == TRUE){
    if (chrom[i,3] == "Tfas") {
      circos.text(chrom$size/2 - n, 2, str_split(name, "_")[[1]][2], sector.index=name,
                  col="olivedrab",cex=0.6, facing = "inside", niceFacing = T)
      n = n + 220000
      if (i > 21 & i < 26){
        n = n + 600000
      }
    } else {
      circos.text(chrom$size - n, 2, str_split(name, "_")[[1]][2], sector.index=name,
                  col="darkgreen",cex=0.6, facing = "inside", niceFacing = T)
      n = n + 20000
      if (i > 38){
        n = n + 2400000
        } 
    }
    } else {
    if ((i > 22 & i < 26) | i > 45){
      if (chrom[i,3] == "Tfas") {
        circos.text(chrom$size/2 - n, 3, str_split(name, "_")[[1]][2], sector.index=name,
                    col="olivedrab",cex=0.6, facing = "inside", niceFacing = T)
        n = n + 600000
      } else {
        circos.text(chrom$size - n, 3, str_split(name, "_")[[1]][2], sector.index=name,
                    col="darkgreen",cex=0.6, facing = "inside", niceFacing = T)
        n = n + 600000
      }
    } else {
      if (chrom[i,3] == "Tfas") {
        circos.text(chrom$size/2 - n , 2, str_split(name, "_")[[1]][2], sector.index=name,
                    col="olivedrab",cex=0.6, facing = "inside", niceFacing = T)
        n = n + 220000
      } else {
        circos.text(chrom$size - n, 2, str_split(name, "_")[[1]][2], sector.index=name,
                    col="darkgreen",cex=0.6, facing = "inside", niceFacing = T)
        n = n + 20000
        if (i > 38){
          n = n + 2400000
        }
      }
    }
  }
}



# Add species lines
draw.sector(91, 290, rou1 = 1.07, rou2 = 1.08, col = "olivedrab", border="NA")
draw.sector(288, 93, rou1 = 1.07, rou2 = 1.08, col = "darkgreen", border="NA")
circos.text(chrom$size/2, 8, "T. fasciculata", sector.index="Tfas_chr11",col="olivedrab",cex=0.9, facing = "bending.inside")
circos.text(chrom$size/2, 8, "T. leiboldiana", sector.index="Tlei_chr11",col="darkgreen",cex=0.9, facing = "bending.inside")

#-------------------TRACK 1: GENE DENSITY-------------------#

## Read in gene content files
# All gene counts
#gene_counts_per_mb_windows <- read.table("Gene_counts_per_1MB_windows.Tfas-Tlei.mainScaffolds.curatedOGs.txt", header = T)
gene_counts_per_mb_windows <- read.table(args[[2]], header = T)
print("Darwing first track: Gene density...")
circos.track(gene_counts_per_mb_windows$chrom, y = gene_counts_per_mb_windows$gene_counts, x = gene_counts_per_mb_windows$start_window, 
              bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "seagreen")
               circos.yaxis(c("left"), sector.index = "Tfas_chr1", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.15, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 50)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

#-------------------TRACK 3: TE DENSITY-------------------#

#TE_content_per_mb_windows <- read.table("TE_content_Tfas-Tlei_per1MB-window_python.txt", header = T, sep = "\t")
TE_content_per_mb_windows <- read.table(args[[3]], header = T,sep = "\t")
print("Drawing second track: TE density...")
circos.track(TE_content_per_mb_windows$chrom, x = TE_content_per_mb_windows$start_window,
             y = TE_content_per_mb_windows$perc_te, bg.col = "grey92", ylim = c(0,100),
             panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "orange")
               circos.yaxis(c("left"), sector.index = "Tfas_chr1", labels = F, at = c(0,10,20),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.17, bg.border = "black")

for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 20)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

#-------------------TRACK 3: DE genes-------------------#

#DE_genes <- read.table("DE_genes_Tfas.txt", header = T)
print("Drawing third track: DE genes...")
DE_genes <- read.table(args[[3]], header = T)
DE_genes <- DE_genes[,c(2:4)]
DE_genes$value <- 1
colnames(DE_genes) <- c("chr", "start", "end", "value1")

circos.genomicTrack(DE_genes, ylim = c(0,10), panel.fun = function(region,value,...)  {
  circos.genomicRect(region, value, ytop = 10, ybottom = 0, col = "yellow", lwd = .01)
}, track.height = 0.17)


#-------------------TRACK 4: SYNTENY-------------------#

### Add links from Tfas to Tlei
#synteny_genes <- read.table("orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.one-to-one.circlize.txt", 
#                            header = T,sep = "\t")
synteny_genes <- read.table(args[[4]], header = T,sep = "\t")
print("Drawing fourth track: Synteny...")
# Create color palette for links
nb.cols <- 25
mycolors <- sample(colorRampPalette(brewer.pal(8, "Set1"))(nb.cols))

### Add links from Tfas to Tlei
for (j in 1:25){
  loc=paste0("Tfas_chr", j)
  genes <- synteny_genes[synteny_genes$Tfas_chrom == loc,]
  for (i in 1:nrow(genes)){
    circos.link(sector.index1=genes[i,3], genes[i,4], sector.index2=genes[i,7], 
                genes[i,8],col=mycolors[j], lwd = .1)
  }
}

dev.off()
print("Done!")