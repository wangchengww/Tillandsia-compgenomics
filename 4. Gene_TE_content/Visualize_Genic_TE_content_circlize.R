## Set working directory and load needed packages
setwd("/Users/clara/bio-info_phd/Comparative_genomics/Circular_figure/")
library(circlize)

## Load a list of chromosomes for both genomes
chrom <- read.table("chromosomes_Tfas_Tlei.txt")
colnames(chrom) <- c("name", "size")
chrom$species <- c(rep(c("Tfas"), 25), rep(c("Tlei"), 26))

## Read in gene content files
# All gene counts
Tfas_all_gene_counts_per_mb_windows <- read.table("Tfas_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
Tlei_all_gene_counts_per_mb_windows <- read.table("Tlei_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
# All gene content (proportional)
Tfas_all_gene_content_per_1_mb_windows <- read.table("Tfas_gene_density_1Mb_windows.all_orthogroups_25_scaffolds.txt")
Tlei_all_gene_content_per_1_mb_windows <- read.table("Tlei_gene_density_1Mb_windows.all_orthogroups_25_scaffolds.txt")
# One-to-one orthologues counts
Tfas_one_gene_counts_per_mb_windows <- read.table("Tfas_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
Tlei_one_gene_counts_per_mb_windows <- read.table("Tlei_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
# One-to-one orthologues content (proportional)
Tfas_one_gene_content_per_1_mb_windows <- read.table("Tfas_gene_density_1Mb_windows.one-to-one_orthogroups_25_scaffolds.txt", header = T)
Tlei_one_gene_content_per_1_mb_windows <- read.table("Tlei_gene_density_1Mb_windows.one-to-one_orthogroups_25_scaffolds.txt", header = T)
# Multi-copy genes counts
Tfas_multi_gene_counts_per_mb_windows <- read.table("Tfas_multicopy_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
Tlei_multi_gene_counts_per_mb_windows <- read.table("Tlei_multicopy_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)

## Read in repetitive content files
Tfas_masked1 <- read.table("TE_content_Tfas_per1MBwindow_py.txt", header = T)
Tfas_masked_percentage_per_mb_windows <- Tfas_masked1
Tfas_masked2 <- read.table("Tfas_TE_density_1Mb_windows.25scaffolds.txt", header = T)
Tfas_masked_percentage_per_mb_windows <- as.data.frame(cbind(Tfas_masked2$chrom, Tfas_masked1$start_window,
                                               Tfas_masked1$end_window, Tfas_masked1$perc_te))
colnames(Tfas_masked_percentage_per_mb_windows) <- c("chrom", "start_window", "end_window", "percentage_te")

Tlei_masked1 <- read.table("TE_content_Tlei_per1MBwindow_py.txt", header = T)
Tlei_masked2 <- read.table("Tlei_TE_density_1Mb_windows.25scaffolds.txt", header = T)
Tlei_masked_percentage_per_mb_windows <- as.data.frame(cbind(Tlei_masked2$chrom, Tlei_masked1$start_window,
                                                             Tlei_masked1$end_window, Tlei_masked1$perc_te))
colnames(Tlei_masked_percentage_per_mb_windows) <- c("chrom", "start_window", "end_window", "percentage_te")


#---------------------------------T.FASCICULATA------------------------------


## Make plot
circos.clear()
# Plot initialization
circos.par("track.height"=0.05, cell.padding=c(0, 0, 0, 0),
           start.degree=96,gap.degree=5,points.overflow.warning=FALSE,
           track.margin=c(0.01,0.01))
pos <- cbind(rep(0, 25), Tfas_chrom$size)
circos.initialize(sectors = Tfas_chrom$name, xlim = pos)
# Add blocks representing the chromosomes
circos.track(Tfas_chrom$name, ylim = c(0, 1))

# Add chromosome names
circos.text(Tfas_chrom$size/2, 2, "chr1", sector.index="Scaffold_2199",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr2", sector.index="Scaffold_658",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr3", sector.index="Scaffold_2318",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr4", sector.index="Scaffold_2041",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr5", sector.index="Scaffold_2096",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr6", sector.index="Scaffold_2315",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr7", sector.index="Scaffold_2317",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr8", sector.index="Scaffold_2287",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr9", sector.index="Scaffold_326",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr10", sector.index="Scaffold_1174",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr11", sector.index="Scaffold_2314",col="olivedrab",cex=0.45, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr12", sector.index="Scaffold_2313",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/2, 2, "chr13", sector.index="Scaffold_2225",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/3, 2, "chr14", sector.index="Scaffold_1383",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr15", sector.index="Scaffold_819",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr16", sector.index="Scaffold_15",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr17", sector.index="Scaffold_2073",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr18", sector.index="Scaffold_2321",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr19", sector.index="Scaffold_2320",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr20", sector.index="Scaffold_346",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr21", sector.index="Scaffold_845",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr22", sector.index="Scaffold_2316",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr23", sector.index="Scaffold_2312",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr24", sector.index="Scaffold_53",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tfas_chrom$size/4, 2, "chr25", sector.index="Scaffold_2319",col="olivedrab",cex=0.4, facing = "inside")

# Add species name
circos.text(Tfas_chrom$size/2, 5, "T. fasciculata", sector.index="Scaffold_2096",col="olivedrab",cex=0.8, facing = "bending.inside")

# Gene density all
circos.track(Tfas_all_gene_content_per_1_mb_windows$chrom, x = as.numeric(Tfas_all_gene_content_per_1_mb_windows$V2),
             y = as.numeric(Tfas_all_gene_content_per_1_mb_windows$V4), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "seagreen")
               circos.yaxis(c("left"), sector.index = "Scaffold_2199", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 50)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

# Gene counts all
circos.track(Tfas_all_gene_counts_per_mb_windows$chrom, x = as.numeric(Tfas_all_gene_counts_per_mb_windows$V2),
             y = as.numeric(Tfas_all_gene_counts_per_mb_windows$counts), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "seagreen")
               circos.yaxis(c("left"), sector.index = "Scaffold_2199", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 50)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}
# Gene density one to one
circos.track(Tfas_one_gene_content_per_1_mb_windows$chrom, x = as.numeric(Tfas_one_gene_content_per_1_mb_windows$V2),
             y = as.numeric(Tfas_one_gene_content_per_1_mb_windows$percentage_genic), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "seagreen")
               circos.yaxis(c("left"), sector.index = "Scaffold_2199", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 50)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

# Gene counts one to one
circos.track(Tfas_one_gene_counts_per_mb_windows$chrom, x = as.numeric(Tfas_one_gene_counts_per_mb_windows$start_window),
             y = as.numeric(Tfas_one_gene_counts_per_mb_windows$gene_count), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "palegreen3")
               circos.yaxis(c("left"), sector.index = "Scaffold_2199", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 20)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

# Multicopy gene counts
circos.track(Tfas_multi_gene_counts_per_mb_windows$chrom, x = as.numeric(Tfas_multi_gene_counts_per_mb_windows$start_window),
             y = as.numeric(Tfas_multi_gene_counts_per_mb_windows$multicopy_gene_count), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "tomato2")
               circos.yaxis(c("left"), sector.index = "Scaffold_2199", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 20)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

# TE density
circos.track(Tfas_masked_percentage_per_mb_windows$chrom, x = as.numeric(Tfas_masked_percentage_per_mb_windows$start_window),
             y = as.numeric(Tfas_masked_percentage_per_mb_windows$percentage_te), bg.col = "grey92", ylim = c(0,100),
             panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "orange")
               circos.yaxis(c("left"), sector.index = "Scaffold_2199", labels = F, at = c(0,10,20),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")

for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 20)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}


#---------------------------------T.LEIBOLDIANA------------------------------

# Make plot
circos.clear()
# Plot initialization
circos.par("track.height"=0.05, cell.padding=c(0, 0, 0, 0),
           start.degree=96,gap.degree=5,points.overflow.warning=FALSE,
           track.margin=c(0.01,0.01))
pos <- cbind(rep(0, 25), Tlei_chrom$size)
circos.initialize(sectors = Tlei_chrom$name, xlim = pos)
# Add blocks representing the chromosomes
circos.track(Tlei_chrom$name, ylim = c(0, 1))

# Add chromosome names
circos.text(Tlei_chrom$size/2, 2, "chr1", sector.index="Scaffold_8399",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr2", sector.index="Scaffold_10370",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr3", sector.index="Scaffold_247",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr4", sector.index="Scaffold_4655",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr5", sector.index="Scaffold_10302",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr6", sector.index="Scaffold_10429",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr7", sector.index="Scaffold_5555",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr8", sector.index="Scaffold_9938",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr9", sector.index="Scaffold_4582",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr10", sector.index="Scaffold_7805",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr11", sector.index="Scaffold_269",col="olivedrab",cex=0.45, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr12", sector.index="Scaffold_6112",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/2, 2, "chr13", sector.index="Scaffold_10423",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/3, 2, "chr14", sector.index="Scaffold_7317",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/4, 2, "chr15", sector.index="Scaffold_4198",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/4, 2, "chr16", sector.index="Scaffold_7315",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/4, 2, "chr17", sector.index="Scaffold_10425",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/4, 2, "chr18", sector.index="Scaffold_10428",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/4, 2, "chr19", sector.index="Scaffold_10433",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/5, 2, "chr20", sector.index="Scaffold_612",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/5, 2, "chr21", sector.index="Scaffold_10430",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/5, 2, "chr22", sector.index="Scaffold_10431",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/5, 2, "chr23", sector.index="Scaffold_10426",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/5, 2, "chr24", sector.index="Scaffold_10432",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/5, 2, "chr25", sector.index="Scaffold_10424",col="olivedrab",cex=0.4, facing = "inside")
circos.text(Tlei_chrom$size/5, 2, "chr26", sector.index="Scaffold_10427",col="olivedrab",cex=0.4, facing = "inside")

# Add species name
circos.text(Tlei_chrom$size/2, 5, "T. leiboldiana", sector.index="Scaffold_10302",col="olivedrab",cex=0.8, facing = "bending.inside")

# Gene density all
circos.track(Tlei_all_gene_content_per_1_mb_windows$chrom, x = as.numeric(Tlei_all_gene_content_per_1_mb_windows$V2),
             y = as.numeric(Tlei_all_gene_content_per_1_mb_windows$V4), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "seagreen")
               circos.yaxis(c("left"), sector.index = "Scaffold_8399", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 50)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

# Gene counts all
circos.track(Tlei_all_gene_counts_per_mb_windows$chrom, x = as.numeric(Tlei_all_gene_counts_per_mb_windows$V2),
             y = as.numeric(Tlei_all_gene_counts_per_mb_windows$counts), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "seagreen")
               circos.yaxis(c("left"), sector.index = "Scaffold_8399", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 50)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}
# Gene density one to one
circos.track(Tlei_one_gene_content_per_1_mb_windows$chrom, x = as.numeric(Tlei_one_gene_content_per_1_mb_windows$V2),
             y = as.numeric(Tlei_one_gene_content_per_1_mb_windows$percentage_genic), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "seagreen")
               circos.yaxis(c("left"), sector.index = "Scaffold_8399", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 50)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

# Gene counts one to one
circos.track(Tlei_one_gene_counts_per_mb_windows$chrom, x = as.numeric(Tlei_one_gene_counts_per_mb_windows$start_window),
             y = as.numeric(Tlei_one_gene_counts_per_mb_windows$gene_count), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "palegreen3")
               circos.yaxis(c("left"), sector.index = "Scaffold_8399", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 20)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

# Multicopy gene counts
circos.track(Tlei_multi_gene_counts_per_mb_windows$chrom, x = as.numeric(Tlei_multi_gene_counts_per_mb_windows$start_window),
             y = as.numeric(Tlei_multi_gene_counts_per_mb_windows$multicopy_gene_count), bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "tomato2")
               circos.yaxis(c("left"), sector.index = "Scaffold_8399", labels = F, at = c(0,70,146),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 20)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

# TE density
circos.track(Tlei_masked_percentage_per_mb_windows$chrom, x = as.numeric(Tlei_masked_percentage_per_mb_windows$start_window),
             y = as.numeric(Tlei_masked_percentage_per_mb_windows$percentage_te), bg.col = "grey92", ylim = c(0,100),
             panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = "orange")
               circos.yaxis(c("left"), sector.index = "Scaffold_8399", labels = F, at = c(0,10,20),
                            labels.cex = 0.3, labels.col="khaki4", tick.length = 2)
             }, track.height = 0.2, bg.border = "black")

for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 20)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}
