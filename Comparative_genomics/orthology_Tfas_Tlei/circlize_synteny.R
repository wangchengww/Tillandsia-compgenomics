## Set working directory and load needed packages
setwd("/Users/clara/bio-info_phd/Comparative_genomics/orthology_Tfas_Tlei/synteny_25_scaffolds/")
library(circlize)

## Load a list of chromosomes for both genomes
chrom <- read.table("chromosomes_Tfas_Tlei.txt")
colnames(chrom) <- c("name", "size")
chrom$species <- c(rep(c("Tfas"), 25), rep(c("Tlei"), 26))
#Create list of start and end positions of each chromosome
start <- c()
end <- c()
b = 0
for (i in 1:nrow(chrom)){
  if (b == 0) {
    s <- 0
    e <- chrom[i,2]
    start <- c(start, s)
    end <- c(end, e)
    b = b +1
  } else {
    s <- e + 1
    e <- e + chrom[i,2]
    start <- c(start, s)
    end <- c(end, e)
  }
}
chrom$start <- start
chrom$end <- end
# Make matrix of start and end position to initialize the circular plot
pos <- cbind(rep(0, 51), chrom$size)

# Create color palette for links
library(RColorBrewer)
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

# Make plot
circos.clear()
# Plot initialization
circos.par("track.height"=0.05, cell.padding=c(0.02, 0, 0.02, 0), 
           start.degree=90,gap.degree=4,points.overflow.warning=FALSE, 
           track.margin=c(0,0))
circos.initialize(sectors = chrom$name, xlim = pos) 

# Add blocks representing the chromosomes
circos.track(chrom$name, ylim = c(0, 1))

# Add chromosome names
circos.text(chrom$size/2, 3, "chr1", sector.index="Scaffold_2199",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr2", sector.index="Scaffold_658",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr3", sector.index="Scaffold_2318",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr4", sector.index="Scaffold_2041",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr5", sector.index="Scaffold_2096",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr6", sector.index="Scaffold_2315",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr7", sector.index="Scaffold_2317",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr8", sector.index="Scaffold_2287",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr9", sector.index="Scaffold_326",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr10", sector.index="Scaffold_1174",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr11", sector.index="Scaffold_2314",col="olivedrab",cex=0.45, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr12", sector.index="Scaffold_2313",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr13", sector.index="Scaffold_2225",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr14", sector.index="Scaffold_1383",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr15", sector.index="Scaffold_819",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr16", sector.index="Scaffold_15",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr17", sector.index="Scaffold_2073",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr18", sector.index="Scaffold_2321",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr19", sector.index="Scaffold_2320",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr20", sector.index="Scaffold_346",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr21", sector.index="Scaffold_845",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr22", sector.index="Scaffold_2316",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr23", sector.index="Scaffold_2312",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr24", sector.index="Scaffold_53",col="olivedrab",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr25", sector.index="Scaffold_2319",col="olivedrab",cex=0.4, facing = "inside")

circos.text(chrom$size/2, 2.5, "chr1", sector.index="Scaffold_8399",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr2", sector.index="Scaffold_10370",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr3", sector.index="Scaffold_247",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr4", sector.index="Scaffold_4655",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr5", sector.index="Scaffold_10302",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr6", sector.index="Scaffold_10429",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr7", sector.index="Scaffold_5555",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr8", sector.index="Scaffold_9938",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr9", sector.index="Scaffold_4582",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr10", sector.index="Scaffold_7805",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr11", sector.index="Scaffold_269",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr12", sector.index="Scaffold_6112",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr13", sector.index="Scaffold_10423",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr14", sector.index="Scaffold_7317",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr15", sector.index="Scaffold_4198",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr16", sector.index="Scaffold_7315",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr17", sector.index="Scaffold_10425",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr18", sector.index="Scaffold_10428",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr19", sector.index="Scaffold_10433",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr20", sector.index="Scaffold_612",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr21", sector.index="Scaffold_10430",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr22", sector.index="Scaffold_10431",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr23", sector.index="Scaffold_10426",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr24", sector.index="Scaffold_10432",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 2.5, "chr25", sector.index="Scaffold_10424",col="darkgreen",cex=0.4, facing = "inside")
circos.text(chrom$size/2, 3, "chr26", sector.index="Scaffold_10427",col="darkgreen",cex=0.4, facing = "inside")

# Add species lines
draw.sector(91, 290, rou1 = 1.07, rou2 = 1.08, col = "olivedrab", border="NA")
draw.sector(288, 93, rou1 = 1.07, rou2 = 1.08, col = "darkgreen", border="NA")
circos.text(chrom$size/2, 8, "T. fasciculata", sector.index="Scaffold_2314",col="olivedrab",cex=0.8, facing = "bending.inside")
circos.text(chrom$size/2, 8, "T. leiboldiana", sector.index="Scaffold_269",col="darkgreen",cex=0.8, facing = "bending.inside")

### Add links from Tfas to Tlei
homeo=read.table("chr1_2199_Tfas_25chrom.txt",h=F, sep="\t")
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[1], lwd = .1)
}
homeo=read.table("chr2_658_Tfas_25chrom.txt",h=F, sep='\t')
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[13], lwd = .1)
}
homeo=read.table("chr3_2318_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[2], lwd = .1)
}
homeo=read.table("chr4_2041_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[14], lwd = .1)
}
homeo=read.table("chr5_2096_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[3], lwd = .1)
}
homeo=read.table("chr6_2315_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[15], lwd = .1)
}
homeo=read.table("chr7_2317_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[4], lwd = .1)
}
homeo=read.table("chr8_2287_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[16], lwd = .1)
}
homeo=read.table("chr9_326_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[5], lwd = .1)
}
homeo=read.table("chr10_1174_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[17], lwd = .1)
}
homeo=read.table("chr11_2314_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[6], lwd = .1)
}
homeo=read.table("chr12_2313_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[18], lwd = .1)
}
homeo=read.table("chr13_2225_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[7], lwd = .1)
}
homeo=read.table("chr14_1383_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[19], lwd = .1)
}
homeo=read.table("chr15_819_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[8], lwd = .1)
}
homeo=read.table("chr16_15_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[20], lwd = .1)
}
homeo=read.table("chr17_2073_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[9], lwd = .1)
}
homeo=read.table("chr18_2321_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[21], lwd = .1)
}
homeo=read.table("chr19_2320_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[10], lwd = .1)
}
homeo=read.table("chr20_346_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[22], lwd = .1)
}
homeo=read.table("chr21_845_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[11], lwd = .1)
}
homeo=read.table("chr22_2316_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[23], lwd = .1)
}
homeo=read.table("chr23_2312_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[12], lwd = .1)
}
homeo=read.table("chr24_53_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[24], lwd = .1)
}
homeo=read.table("chr25_2319_Tfas_25chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,2], homeo[i,3], sector.index2=homeo[i,6], homeo[i,7],col=mycolors[25], lwd = .1)
}

### Color palette for Tlei
nb.cols <- 26
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

### Add links from Tlei to Tfas
homeo=read.table("chr1_8399_Tlei_26chrom.txt",h=F, sep="\t")
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[1], lwd = .5)
}
homeo=read.table("chr2_10370_Tlei_26chrom.txt",h=F, sep='\t')
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[14], lwd = .5)
}
homeo=read.table("chr3_247_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[2], lwd = .5)
}
homeo=read.table("chr4_4655_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[15], lwd = .5)
}
homeo=read.table("chr5_10302_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[3], lwd = .5)
}
homeo=read.table("chr6_10429_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[16], lwd = .5)
}
homeo=read.table("chr7_5555_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[4], lwd = .5)
}
homeo=read.table("chr8_9938_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[17], lwd = .5)
}
homeo=read.table("chr9_4582_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[5], lwd = .5)
}
homeo=read.table("chr10_7805_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[18], lwd = .5)
}
homeo=read.table("chr11_269_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[6], lwd = .5)
}
homeo=read.table("chr12_6112_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[19], lwd = .5)
}
homeo=read.table("chr13_10423_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[7], lwd = .5)
}
homeo=read.table("chr14_7317_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[20], lwd = .5)
}
homeo=read.table("chr15_4198_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[8], lwd = .5)
}
homeo=read.table("chr16_7315_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[21], lwd = .5)
}
homeo=read.table("chr17_10425_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[9], lwd = .5)
}
homeo=read.table("chr18_10428_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[22], lwd = .5)
}
homeo=read.table("chr19_10433_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[10], lwd = .5)
}
homeo=read.table("chr20_612_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[23], lwd = .5)
}
homeo=read.table("chr21_10430_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[11], lwd = .5)
}
homeo=read.table("chr22_10431_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[24], lwd = .5)
}
homeo=read.table("chr23_10426_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[12], lwd = .5)
}
homeo=read.table("chr24_10432_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[25], lwd = .5)
}
homeo=read.table("chr25_10424_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[13], lwd = .5)
}
homeo=read.table("chr26_10427_Tlei_26chrom.txt",h=F)
for (i in 1:nrow(homeo)){
  circos.link(sector.index1=homeo[i,6], homeo[i,7], sector.index2=homeo[i,2], homeo[i,3],col=mycolors[26], lwd = .5)
}

