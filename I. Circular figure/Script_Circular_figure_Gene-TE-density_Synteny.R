#!/usr/bin/Rscript --vanilla

# Code for large circular figure - Genome Paper T.fasciculata / T. leiboldiana
# This Rcode creates a circular figure with the 25 / 26 main scaffolds of both assemblies as halves of the circle. The circle has 3 tracks, a Gene density, a TE density and a synteny plot

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("circlize", "stringr", "RColorBrewer")

#setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure")
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure')

# Load arguments
# 1 is chromosome list, 2 is the gene density, 3 is tTE density, 4 is DE genes, 5 is synteny,
# 6 is output name
args <- commandArgs(trailingOnly = TRUE)
output_name <- args[[6]]
# Read in the complete list of chromosomes with all needed info:
chrom <- read.table("chromosomes_Tfas_Tlei.coordinates-circle.txt.no-Tlei_chr2526", header = T, sep = "\t")
chrom <- read.table(args[[1]], header = T, sep = "\t")

# Make matrix of start and end position to initialize the circular plot
pos <- cbind(rep(0, dim(chrom)[1]), chrom$size)

# Make plot
pdf(paste0(output_name, ".pdf"), width = 10, height = 8)
print("Initializing the plot...")
circos.clear()
# Plot initialization
circos.par("track.height"=0.05, cell.padding=c(0.02, 0, 0.02, 0),
           start.degree=90,gap.degree=1.5,points.overflow.warning=FALSE,
           track.margin=c(0,0))
circos.initialize(sectors = chrom$name, xlim = pos)


# Add blocks representing the chromosomes
circos.track(chrom$name, ylim = c(0, 1), bg.border = "white", track.height = .05)

# Add chromosome names
is.even <- function(x) x %% 2 == 0
n = 0
m = 0
for (i in 1:nrow(chrom)){
  name=chrom[i,1]
  if (is.even(i) == TRUE){
    if (chrom[i,3] == "Tfas") {
      circos.text(chrom$size/2 - n, 1.1, str_split(name, "_")[[1]][2], sector.index=name,
                  col="olivedrab",cex=0.6, facing = "inside", niceFacing = T)
      n = n + 220000
      if (i > 21 & i < 26){
        n = n + 600000
      }
    } else {
      circos.text(chrom$size - n, 1.1, str_split(name, "_")[[1]][2], sector.index=name,
                  col="darkgreen",cex=0.6, facing = "inside", niceFacing = T)
      n = n + 20000
      if (i > 38){
        n = n + 2400000
      } 
    }
  } else {
    if ((i > 22 & i < 26) | i > 45){
      if (chrom[i,3] == "Tfas") {
        circos.text(chrom$size/2 - n, 2.2, str_split(name, "_")[[1]][2], sector.index=name,
                    col="olivedrab",cex=0.6, facing = "inside", niceFacing = T)
        n = n + 600000
      } else {
        circos.text(chrom$size - n, 2.2, str_split(name, "_")[[1]][2], sector.index=name,
                    col="darkgreen",cex=0.6, facing = "inside", niceFacing = T)
        n = n + 600000
      }
    } else {
      if (chrom[i,3] == "Tfas") {
        circos.text(chrom$size/2 - n , 1.1, str_split(name, "_")[[1]][2], sector.index=name,
                    col="olivedrab",cex=0.6, facing = "inside", niceFacing = T)
        n = n + 220000
      } else {
        circos.text(chrom$size - n, 1.1, str_split(name, "_")[[1]][2], sector.index=name,
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
draw.sector(90, 299, rou1 = 1.05, rou2 = 1.06, col = "olivedrab", border="NA")
draw.sector(298, 91, rou1 = 1.05, rou2 = 1.06, col = "darkgreen", border="NA")
circos.text(chrom$size/2, 6, "T. fasciculata", sector.index="Tfas_chr11",col="olivedrab",cex=0.9, facing = "bending.inside")
circos.text(chrom$size/2, 6, "T. leiboldiana", sector.index="Tlei_chr11",col="darkgreen",cex=0.9, facing = "bending.inside")

#-------------------TRACK 1: GENE DENSITY-------------------#

## Read in gene content files
# All gene counts
gene_counts_per_mb_windows <- read.table("Gene_counts_per_1MB_windows.Tfas-Tlei.mainScaffolds.curatedOGs.txt.no-Tlei_chr2526", header = T)
gene_counts_per_mb_windows <- read.table(args[[2]], header = T)
print("Darwing first track: Gene density...")
color_gene_density <- c(Tfas_chr1 = "palegreen3", Tfas_chr2 = "palegreen3", Tfas_chr3 = "palegreen3",
                        Tfas_chr4 = "palegreen3", Tfas_chr5 = "palegreen3", Tfas_chr6 = "palegreen3",
                        Tfas_chr7 = "palegreen3", Tfas_chr8 = "palegreen3", Tfas_chr9 = "palegreen3",
                        Tfas_chr10 = "palegreen3", Tfas_chr11 = "palegreen3", Tfas_chr12 = "palegreen3",
                        Tfas_chr13 = "palegreen3", Tfas_chr14 = "palegreen3", Tfas_chr15 = "palegreen3",
                        Tfas_chr16 = "palegreen3", Tfas_chr17 = "palegreen3", Tfas_chr18 = "palegreen3",
                        Tfas_chr19 = "palegreen3", Tfas_chr20 = "palegreen3", Tfas_chr21 = "palegreen3",
                        Tfas_chr22 = "palegreen3", Tfas_chr23 = "palegreen3", Tfas_chr24 = "palegreen3",
                        Tfas_chr25 = "palegreen3", Tlei_chr1 = "seagreen", Tlei_chr2 = "seagreen",
                        Tlei_chr3 = "seagreen", Tlei_chr4 = "seagreen", Tlei_chr5 = "seagreen",
                        Tlei_chr6 = "seagreen", Tlei_chr7 = "seagreen", Tlei_chr8 = "seagreen",
                        Tlei_chr9 = "seagreen", Tlei_chr10 = "seagreen", Tlei_chr11 = "seagreen",
                        Tlei_chr12 = "seagreen", Tlei_chr13 = "seagreen", Tlei_chr14 = "seagreen",
                        Tlei_chr15 = "seagreen",Tlei_chr16 = "seagreen", Tlei_chr17 = "seagreen",
                        Tlei_chr18 = "seagreen",Tlei_chr19 = "seagreen", Tlei_chr20 = "seagreen",
                        Tlei_chr21 = "seagreen",Tlei_chr22 = "seagreen", Tlei_chr23 = "seagreen",
                        Tlei_chr24 = "seagreen")
circos.track(gene_counts_per_mb_windows$chrom, y = gene_counts_per_mb_windows$gene_counts,
             x = gene_counts_per_mb_windows$start_window,
             bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = color_gene_density[CELL_META$sector.index])
               circos.yaxis(c("left"), sector.index = "Tfas_chr1", labels = T, at = seq(0, CELL_META$ylim[2], by = 50),
                            labels.cex = 0.15, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.15, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 50)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}
#-------------------TRACK 3: TE DENSITY-------------------#

#TE_content_per_mb_windows <- read.table("TE_content_Tfas-Tlei_per1MB-window_python.txt.no-Tlei_chr2526", header = T, sep = "\t")
TE_content_per_mb_windows <- read.table(args[[3]], header = T,sep = "\t")
print("Drawing second track: TE density...")

color_te_density <- c(Tfas_chr1 = "goldenrod1", Tfas_chr2 = "goldenrod1", Tfas_chr3 = "goldenrod1",
                        Tfas_chr4 = "goldenrod1", Tfas_chr5 = "goldenrod1", Tfas_chr6 = "goldenrod1",
                        Tfas_chr7 = "goldenrod1", Tfas_chr8 = "goldenrod1", Tfas_chr9 = "goldenrod1",
                        Tfas_chr10 = "goldenrod1", Tfas_chr11 = "goldenrod1", Tfas_chr12 = "goldenrod1",
                        Tfas_chr13 = "goldenrod1", Tfas_chr14 = "goldenrod1", Tfas_chr15 = "goldenrod1",
                        Tfas_chr16 = "goldenrod1", Tfas_chr17 = "goldenrod1", Tfas_chr18 = "goldenrod1",
                        Tfas_chr19 = "goldenrod1", Tfas_chr20 = "goldenrod1", Tfas_chr21 = "goldenrod1",
                        Tfas_chr22 = "goldenrod1", Tfas_chr23 = "goldenrod1", Tfas_chr24 = "goldenrod1",
                        Tfas_chr25 = "goldenrod1", Tlei_chr1 = "orange2", Tlei_chr2 = "orange2",
                        Tlei_chr3 = "orange2", Tlei_chr4 = "orange2", Tlei_chr5 = "orange2",
                        Tlei_chr6 = "orange2", Tlei_chr7 = "orange2", Tlei_chr8 = "orange2",
                        Tlei_chr9 = "orange2", Tlei_chr10 = "orange2", Tlei_chr11 = "orange2",
                        Tlei_chr12 = "orange2", Tlei_chr13 = "orange2", Tlei_chr14 = "orange2",
                        Tlei_chr15 = "orange2",Tlei_chr16 = "orange2", Tlei_chr17 = "orange2",
                        Tlei_chr18 = "orange2",Tlei_chr19 = "orange2", Tlei_chr20 = "orange2",
                        Tlei_chr21 = "orange2",Tlei_chr22 = "orange2", Tlei_chr23 = "orange2",
                        Tlei_chr24 = "orange2")

circos.track(TE_content_per_mb_windows$chrom, x = TE_content_per_mb_windows$start_window,
             y = TE_content_per_mb_windows$perc_te, bg.col = "grey92", ylim = c(0,100),
             panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = color_te_density[CELL_META$sector.index])
               circos.yaxis(c("left"), sector.index = "Tfas_chr1", labels = T, at = seq(0, CELL_META$ylim[2], by = 20),
                            labels.cex = 0.15, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.17, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 20)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

#-------------------TRACK 3: DE genes-------------------#

#DE_genes <- read.table("DE_genes_Tfas-Tlei.txt.no-Tlei_chr2526", header = T)
print("Drawing third track: DE genes...")
DE_genes <- read.table(args[[4]], header = T)
DE_genes <- DE_genes[,c(2:4)]

n_tfas_genes = sum(grepl("Tfas", DE_genes$chr))
n_tlei_genes = sum(grepl("Tlei", DE_genes$chr))
DE_genes$value <- c(rep(1,n_tfas_genes),rep(2,n_tlei_genes))
colnames(DE_genes) <- c("chr", "start", "end", "value1")
circos.genomicTrack(DE_genes, ylim = c(0,10), panel.fun = function(region,value,...)  {
  i = getI(...)
  circos.genomicRect(region, value, ytop = 10, ybottom = 0, col = ifelse(value[[1]] == 1, "olivedrab", "darkgreen"), border = ifelse(value[[1]] == 1, "olivedrab", "darkgreen"), lwd = .0000001)
}, track.height = 0.1)


#-------------------TRACK 4: SYNTENY-------------------#

### Add links from Tfas to Tlei
#synteny_genes <- read.table("orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.one-to-one.circlize.txt.no-Tlei_chr2526",
#                            header = T,sep = "\t")
synteny_genes <- read.table(args[[5]], header = T,sep = "\t")
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
                genes[i,8],col=mycolors[j], lwd = .03)
  }
}

dev.off()
print("Done!")
