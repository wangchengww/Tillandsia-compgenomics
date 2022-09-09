# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("circlize", "stringr", "RColorBrewer")

#setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure")
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure/Supp_DE_density/')

# Load arguments
# 1 is chromosome list, 2 is the gene density, 3 is tTE density, 4 is DE genes, 5 is synteny,
# 6 is output name
args <- commandArgs(trailingOnly = TRUE)
output_name <- args[[5]]
# Read in the complete list of chromosomes with all needed info:
chrom <- read.table("../chromosomes_Tfas_Tlei.coordinates-circle.txt", header = T, sep = "\t")
chrom <- read.table(args[[1]], header = T, sep = "\t")
chrom$species <- c(rep(c("Tfas"), 25), rep(c("Tlei"), 26))
Tfas_chrom <- chrom[chrom$species=="Tfas",]
Tlei_chrom <- chrom[chrom$species=="Tlei",]

# Make matrix of start and end position to initialize the circular plot
pos <- cbind(rep(0, dim(Tfas_chrom)[1]), Tfas_chrom$size)

# Make plot
pdf(paste0(output_name, ".pdf"), width = 10, height = 8)
print("Initializing the plot...")
circos.clear()
# Plot initialization
circos.par("track.height"=0.05, cell.padding=c(0.02, 0, 0.02, 0),
           start.degree=90,gap.degree=1.5,points.overflow.warning=FALSE,
           track.margin=c(0,0))
circos.initialize(sectors = Tfas_chrom$name, xlim = pos)


# Add blocks representing the chromosomes
circos.track(Tfas_chrom$name, ylim = c(0, 1), bg.border = "white", track.height = .05)

# Add chromosome names
is.even <- function(x) x %% 2 == 0
n = 0
m = 0
for (i in 1:nrow(Tfas_chrom)){
  name=Tfas_chrom[i,1]
  if (is.even(i) == TRUE){
      circos.text(chrom$size/2 - n, 1.1, str_split(name, "_")[[1]][2], sector.index=name,
                  col="olivedrab",cex=0.7, facing = "inside", niceFacing = T)
      n = n + 220000
      if (i > 21 & i < 26){
        n = n + 600000
      }
  } else {
    if ((i > 23 & i < 26) | i > 45){
      circos.text(chrom$size/2 - n, 2.2, str_split(name, "_")[[1]][2], sector.index=name,
                    col="olivedrab",cex=0.7, facing = "inside", niceFacing = T)
      n = n + 600000
    } else {
      circos.text(chrom$size/2 - n , 1.1, str_split(name, "_")[[1]][2], sector.index=name,
                    col="olivedrab",cex=0.7, facing = "inside", niceFacing = T)
      n = n + 220000
    }
  }
}



# Add species lines
draw.sector(90, 299, rou1 = 1.05, rou2 = 1.06, col = "olivedrab", border="NA")
draw.sector(298, 91, rou1 = 1.05, rou2 = 1.06, col = "darkgreen", border="NA")
circos.text(chrom$size/2, 6, "T. fasciculata", sector.index="Tfas_chr11",col="olivedrab",cex=1, facing = "bending.inside")
circos.text(chrom$size/2, 6, "T. leiboldiana", sector.index="Tlei_chr11",col="darkgreen",cex=1, facing = "bending.inside")

DE_genes <- read.table("DE_genes.mapped_to_Tfas.ROBUSTonly.txt", header = T)
DE_genes_density <- read.table("DE_genes_Tfas-Tlei.txt.no-Tlei_chr2526.PER-WINDOW.txt", header = T, sep = "\t")
#DE_genes_density <- read.table(args[[4]], header = T)
DE_genes <- DE_genes[,c(2:4)]

# Calculate density
#DE_genes_density <- data.frame()
#for (i in 1:nrow(chrom)){
#  chr = chrom[i,1]
#  print(chr)
#  starts_win <- seq(1, chrom[i,2], by = 1000000)
#  n <- length(starts_win)
#  genes = DE_genes[DE_genes$chrom == chr,]
#  for (n in 1:n){
#    seq = seq(starts_win[n], starts_win[n]+999999)
#    counts = 0
#    for (i in 1:nrow(genes)){
#      in_window <- genes[i,3] %in% seq
#      if (in_window == "TRUE"){
#        counts = counts + 1
#      } else {
#        counts = counts
#      }
#    }
#    window_entry <- cbind(chr, starts_win[n], starts_win[n]+999999, counts)
#    DE_genes_density <- rbind(DE_genes_density, window_entry)
#  }
#}
#colnames(DE_genes_density) <- c("chrom", "start_window", "end_window", "DE_counts")
#DE_genes_density$DE_counts <- as.numeric(DE_genes_density$DE_counts)
#DE_genes_density$start_window <- as.numeric(DE_genes_density$start_window)

#DE_genes_density$proportion <- DE_genes_density$DE_counts/gene_counts_per_mb_windows$gene_counts
#DE_genes_density[is.na(DE_genes_density)] <- 0
#write.table(DE_genes_density, file = "DE_genes_Tfas-Tlei.txt.no-Tlei_chr2526.PER-WINDOW.txt", quote = F, sep = "\t")

color_de_density <- c(Tfas_chr1 = "olivedrab3", Tfas_chr2 = "olivedrab3", Tfas_chr3 = "olivedrab3",
                        Tfas_chr4 = "olivedrab3", Tfas_chr5 = "olivedrab3", Tfas_chr6 = "olivedrab3",
                        Tfas_chr7 = "olivedrab3", Tfas_chr8 = "olivedrab3", Tfas_chr9 = "olivedrab3",
                        Tfas_chr10 = "olivedrab3", Tfas_chr11 = "olivedrab3", Tfas_chr12 = "olivedrab3",
                        Tfas_chr13 = "olivedrab3", Tfas_chr14 = "olivedrab3", Tfas_chr15 = "olivedrab3",
                        Tfas_chr16 = "olivedrab3", Tfas_chr17 = "olivedrab3", Tfas_chr18 = "olivedrab3",
                        Tfas_chr19 = "olivedrab3", Tfas_chr20 = "olivedrab3", Tfas_chr21 = "olivedrab3",
                        Tfas_chr22 = "olivedrab3", Tfas_chr23 = "olivedrab3", Tfas_chr24 = "olivedrab3",
                        Tfas_chr25 = "olivedrab3", Tlei_chr1 = "forestgreen", Tlei_chr2 = "forestgreen",
                        Tlei_chr3 = "forestgreen", Tlei_chr4 = "forestgreen", Tlei_chr5 = "forestgreen",
                        Tlei_chr6 = "forestgreen", Tlei_chr7 = "forestgreen", Tlei_chr8 = "forestgreen",
                        Tlei_chr9 = "forestgreen", Tlei_chr10 = "forestgreen", Tlei_chr11 = "forestgreen",
                        Tlei_chr12 = "forestgreen", Tlei_chr13 = "forestgreen", Tlei_chr14 = "forestgreen",
                        Tlei_chr15 = "forestgreen",Tlei_chr16 = "forestgreen", Tlei_chr17 = "forestgreen",
                        Tlei_chr18 = "forestgreen",Tlei_chr19 = "forestgreen", Tlei_chr20 = "forestgreen",
                        Tlei_chr21 = "forestgreen",Tlei_chr22 = "forestgreen", Tlei_chr23 = "forestgreen",
                        Tlei_chr24 = "forestgreen")

circos.track(DE_genes_density$chrom, y = DE_genes_density$proportion,
             x = DE_genes_density$start_window, ylim = c(0,1),
             bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = color_de_density[CELL_META$sector.index],
                            lwd = .8)
               circos.yaxis(c("left"), sector.index = "Tfas_chr1", labels = T, at = seq(0, CELL_META$ylim[2], by = 2.5),
                            labels.cex = 0.25, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.12, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 2.5)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
#}