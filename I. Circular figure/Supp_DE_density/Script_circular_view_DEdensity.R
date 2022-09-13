# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("circlize", "stringr", "RColorBrewer", "gridGraphics", "gridExtra", "cowplot")

#setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure")
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure/Supp_DE_density/')

chrom <- read.table("../chromosomes_Tfas_Tlei.coordinates-circle.txt", header = T, sep = "\t")
chrom$species <- c(rep(c("Tfas"), 25), rep(c("Tlei"), 26))
Tfas_chrom <- chrom[chrom$species=="Tfas",]
Tlei_chrom <- chrom[chrom$species=="Tlei",]

############## CALCULATE DENSITY OF DE GENES ##############

#all_genes <- read.table("../Gene_counts_per_1MB_windows.Tfas-Tlei.mainScaffolds.curatedOGs.txt", header = T)
#all_genes$gene_counts <- as.numeric(all_genes$gene_counts)
#all_genes_Tfas <- all_genes[startsWith(all_genes$chrom, "Tfas"),]
#all_genes_Tlei <- all_genes[startsWith(all_genes$chrom, "Tlei"),]
#DE_genes_Tfas <- read.table("DE_genes.mapped_to_Tfas.ROBUSTonly.txt", header = T)
#DE_genes_Tlei <- read.table("DE_genes.mapped_to_Tlei.ROBUSTonly.txt", header = T)
#DE_genes_Tlei <- DE_genes_Tlei[,c(2:4)]
#DE_genes_density <- data.frame()
#for (i in 1:nrow(Tlei_chrom)){
#  chr = Tlei_chrom[i,1]
#  print(chr)
#  starts_win <- seq(1, Tlei_chrom[i,2], by = 1000000)
#  n <- length(starts_win)
#  genes = DE_genes_Tlei[DE_genes_Tlei$chrom == chr,]
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
#DE_genes_density$proportion <- DE_genes_density$DE_counts/all_genes_Tlei$gene_counts
#DE_genes_density[is.na(DE_genes_density)] <- 0
#write.table(DE_genes_density, file = "DE_genes.mapped_to_Tlei.PER-WINDOW.density.txt", quote = F, sep = "\t")

############## CIRCLE TFASCICULATA ##############

circos_plot <- function(density, chrom, col, Name){
  # Make matrix of start and end position to initialize the circular plot
  pos <- cbind(rep(0, dim(chrom)[1]), chrom$size)
  # Make plot
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
      circos.text(chrom$size/2 - n, 1.1, str_split(name, "_")[[1]][2], sector.index=name,
                  col=col,cex=0.7, facing = "inside", niceFacing = T)
      n = n + 220000
      if (i > 21 & i < 26){
        n = n + 600000
      }
    } else {
      if ((i > 23 & i < 26) | i > 45){
        circos.text(chrom$size/2 - n, 2.2, str_split(name, "_")[[1]][2], sector.index=name,
                    col=col,cex=0.7, facing = "inside", niceFacing = T)
        n = n + 600000
      } else {
        circos.text(chrom$size/2 - n , 1.1, str_split(name, "_")[[1]][2], sector.index=name,
                    col=col,cex=0.7, facing = "inside", niceFacing = T)
        n = n + 220000
      }
    }
  }
  # Add species name in the center
  text(0, 0, Name, cex = .8, col = col)
  # Read in DE gene density data
  DE_genes_density <- read.table(density,header = T, sep = "\t")
  # Plot track
  circos.track(DE_genes_density$chrom, y = DE_genes_density$proportion,
               x = DE_genes_density$start_window, ylim = c(-.2,1),
               bg.col = "grey92", panel.fun = function(x, y) {
                 circos.lines(x, y, area = T, col = col, lwd = .8)
                 circos.yaxis(c("left"), sector.index = DE_genes_density$chrom[1], labels = T, at = seq(0, CELL_META$ylim[2], by = 2.5),
                              labels.cex = 0.25, labels.col="khaki4", tick.length = 2)
               }, track.height = 0.12, bg.border = "black")
  for(sn in get.all.sector.index()) {
    set.current.cell(sector.index = sn, track.index = get.current.track.index())
    breaks = seq(0, CELL_META$ylim[2], by = 2.5)
    for(b in breaks) {
      circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
     }
   }
 }

p1 <- function() circos_plot("DE_genes.mapped_to_Tfas.PER-WINDOW.density.txt", 
                             Tfas_chrom, "olivedrab", "T. fasciculata")
p2 <- function() circos_plot("DE_genes.mapped_to_Tlei.PER-WINDOW.density.txt", 
                             Tlei_chrom, "darkgreen", "T. leiboldiana")
cowplot::plot_grid(p2)


############## CIRCLE TLEIBOLDIANA ##############

# Make matrix of start and end position to initialize the circular plot
pos <- cbind(rep(0, dim(Tlei_chrom)[1]), Tlei_chrom$size)
# Make plot
circos.clear()
# Plot initialization
circos.par("track.height"=0.05, cell.padding=c(0.02, 0, 0.02, 0),
           start.degree=90,gap.degree=1.5,points.overflow.warning=FALSE,
           track.margin=c(0,0))
circos.initialize(sectors = Tlei_chrom$name, xlim = pos)

# Add blocks representing the chromosomes
circos.track(Tlei_chrom$name, ylim = c(0, 1), bg.border = "white", track.height = .05)
# Add chromosome names
is.even <- function(x) x %% 2 == 0
n = 0
m = 0
for (i in 1:nrow(Tlei_chrom)){
  name=Tlei_chrom[i,1]
  if (is.even(i) == TRUE){
    circos.text(Tlei_chrom$size/2 - n, 1.1, str_split(name, "_")[[1]][2], sector.index=name,
                col="olivedrab",cex=0.7, facing = "inside", niceFacing = T)
    n = n + 220000
    if (i > 21 & i < 26){
      n = n + 600000
    }
  } else {
    if ((i > 23 & i < 26) | i > 45){
      circos.text(Tlei_chrom$size/2 - n, 2.2, str_split(name, "_")[[1]][2], sector.index=name,
                  col="olivedrab",cex=0.7, facing = "inside", niceFacing = T)
      n = n + 600000
    } else {
      circos.text(Tlei_chrom$size/2 - n , 1.1, str_split(name, "_")[[1]][2], sector.index=name,
                  col="olivedrab",cex=0.7, facing = "inside", niceFacing = T)
      n = n + 220000
    }
  }
}
# Add species name in the center
text(0, 0, "T. lei", cex = .8, col = "olivedrab")
# Read in DE gene density data
DE_genes_density <- read.table("DE_genes.mapped_to_Tlei.PER-WINDOW.density.txt",
                                    header = T, sep = "\t")
# Plot track
circos.track(DE_genes_density$chrom, y = DE_genes_density$proportion,
             x = DE_genes_density$start_window, ylim = c(-.2,1),
             bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "olivedrab", lwd = .8)
               circos.yaxis(c("left"), sector.index = DE_genes_density$chrom[1], labels = T, at = seq(0, CELL_META$ylim[2], by = 2.5),
                            labels.cex = 0.25, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.12, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 2.5)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}
pdf(paste0(output_name, ".pdf"), width = 10, height = 8)
