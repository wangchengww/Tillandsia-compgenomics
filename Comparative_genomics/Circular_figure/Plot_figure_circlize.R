## Set working directory and load needed packages
setwd("/Users/clara/bio-info_phd/Comparative_genomics/Circular_figure/")
library(circlize)

## Load a list of chromosomes for both genomes
chrom <- read.table("chromosomes_Tfas_Tlei copy.txt")
colnames(chrom) <- c("name", "size")
chrom$species <- c(rep(c("Tfas"), 25), rep(c("Tlei"), 26))

## Load in list of curated orthologous genes
genes_table <- read.delim("../orthology_Tfas_Tlei/orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt", sep = "\t", header = F)

### Make circle for each ref genome, starting with Tfas

#1st track: GENE DENSITY - ALL GENES

# Select T.fas genes and remove all species specific orthogroups
Tfas_genes <- subset(genes_table, grepl("^Tfasc.v1",V1))
Tfas_genes <- Tfas_genes[!(Tfas_genes$V8 == 0 & Tfas_genes$V10 == 0),]

# Calculate gene density
Tfas_chrom <- chrom[chrom$species=="Tfas",]

#Tfas_genes <- read.table("Tfas_orthology_info_per_scaffold.run2.txt")
Tfas_genes <- Tfas_genes[,c(1,2,3,4)]
colnames(Tfas_genes) <- c("gene_id", "chrom", "start", "end")

## Loop to create per 1 Mb window count of starting genes
Tfas_all_gene_counts_per_mb_windows <- data.frame()
for (i in 1:nrow(Tfas_chrom)){
  chrom = Tfas_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tfas_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  genes = Tfas_genes[Tfas_genes$chrom == chrom,]
  for (n in 1:n){
    seq = seq(starts_win[n], starts_win[n]+999999)
    counts = 0
    for (i in 1:nrow(genes)){
      in_window <- genes[i,3] %in% seq
      if (in_window == "TRUE"){
        counts = counts + 1
      } else {
        counts = counts
      }
    }
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, counts)
    Tfas_all_gene_counts_per_mb_windows <- rbind(Tfas_all_gene_counts_per_mb_windows, window_entry)
  }
}
#write.table(x = counts_per_mb_windows, file = "Tfas_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")
Tfas_all_gene_counts_per_mb_windows <- read.table("Tfas_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
summary(Tfas_all_gene_counts_per_mb_windows)

## Loop to create proportion of genic bases per 1 Mb window
Tfas_all_gene_content_per_1_mb_windows <- data.frame()
rest <- data.frame()
for (i in 1:nrow(Tfas_chrom)){
  chrom = Tfas_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tfas_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  chrom_genes = Tfas_genes[Tfas_genes$chrom == chrom,]
  for (n in 1:n){
    start <- starts_win[n]
    end <- starts_win[n]+1000000
    window_genes <- chrom_genes[chrom_genes$start > start & chrom_genes$start < end,]
    window_genes$length <- window_genes$end - window_genes$start
    genic_length <- sum(window_genes$length)
    rest_length <- rest$end - start
    if (dim(rest)[1] != 0){
      genic_length <- genic_length + rest_length
    }
    percentage_genic <- (genic_length/1000000)*100
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, as.numeric(percentage_genic))
    Tfas_all_gene_content_per_1_mb_windows <- rbind(Tfas_all_gene_content_per_1_mb_windows, window_entry)
    rest <- window_genes[window_genes$end > end,]
  }
}
#write.table(Tfas_all_gene_content_per_1_mb_windows, file = "Tfas_gene_density_1Mb_windows.all_orthogroups_25_scaffolds.txt")
Tfas_all_gene_content_per_1_mb_windows <- read.table("Tfas_gene_density_1Mb_windows.all_orthogroups_25_scaffolds.txt")

## 2nd track: GENE DENSITY - 1:1 ORTHOGROUPS

# Select T.fas genes with 1:1 relationship
Tfas_genes_one <- Tfas_genes[(Tfas_genes$V9 == 1 & Tfas_genes$V10 == 1),]

# Calculate gene density
Tfas_chrom <- chrom[chrom$species=="Tfas",]

#Tfas_genes <- read.table("Tfas_orthology_info_per_scaffold.run2.txt")
Tfas_genes_one <- Tfas_genes_one[,c(1,2,3,4)]
colnames(Tfas_genes_one) <- c("gene_id", "chrom", "start", "end")

## Loop to create per 1 Mb window count of starting genes
Tfas_one_genecounts_per_mb_windows <- data.frame()
for (i in 1:nrow(Tfas_chrom)){
  chrom = Tfas_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tfas_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  genes = Tfas_genes_one[Tfas_genes_one$chrom == chrom,]
  for (n in 1:n){
    seq = seq(starts_win[n], starts_win[n]+999999)
    counts = 0
    for (i in 1:nrow(genes)){
      in_window <- genes[i,3] %in% seq
      if (in_window == "TRUE"){
        counts = counts + 1
      } else {
        counts = counts
      }
    }
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, as.numeric(counts))
    Tfas_one_gene_counts_per_mb_windows <- rbind(Tfas_one_genecounts_per_mb_windows, window_entry)
  }
}
#write.table(x = Tfas_one_gene_counts_per_mb_windows, file = "Tfas_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")
Tfas_one_gene_counts_per_mb_windows <- read.table("Tfas_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
colnames(Tfas_one_gene_counts_per_mb_windows) <- c("chrom", "start_window", "end_window", "gene_count")
summary(Tfas_one_gene_counts_per_mb_windows)

## Loop to create proportion of genic bases per 1 Mb window
Tfas_one_gene_content_per_1_mb_windows <- data.frame()
rest <- data.frame()
for (i in 1:nrow(Tfas_chrom)){
  chrom = Tfas_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tfas_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  chrom_genes = Tfas_genes_one[Tfas_genes_one$chrom == chrom,]
  for (n in 1:n){
    start <- starts_win[n]
    end <- starts_win[n]+1000000
    window_genes <- chrom_genes[chrom_genes$start > start & chrom_genes$start < end,]
    window_genes$length <- window_genes$end - window_genes$start
    genic_length <- sum(window_genes$length)
    print(genic_length)
    rest_length <- rest$end - start
    if (dim(rest)[1] != 0){
      print("There's a rest! New length:")
      genic_length <- genic_length + rest_length
      print(genic_length)
    }
    percentage_genic <- (genic_length/1000000)*100
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, percentage_genic)
    Tfas_one_gene_content_per_1_mb_windows <- rbind(Tfas_one_gene_content_per_1_mb_windows, window_entry)
    rest <- window_genes[window_genes$end > end,]
  }
}
#write.table(Tfas_one_gene_content_per_1_mb_windows, file = "Tfas_gene_density_1Mb_windows.one-to-one_orthogroups_25_scaffolds.txt")
Tfas_one_gene_content_per_1_mb_windows <- read.table("Tfas_gene_density_1Mb_windows.one-to-one_orthogroups_25_scaffolds.txt", header = T)
summary(Tfas_one_gene_content_per_1_mb_windows)

## 3nd track: MULTICOPY genes in Tfas

Tfas_genes_multi <- Tfas_genes[!(Tfas_genes$V9 == 1),]
Tfas_genes_multi <- Tfas_genes_multi[,c(1,2,3,4)]
colnames(Tfas_genes_multi) <- c("gene_id", "chrom", "start", "end")
## Loop to create per 1 Mb window count of starting genes
Tfas_multi_gene_counts_per_mb_windows <- data.frame()
for (i in 1:nrow(Tfas_chrom)){
  chrom = Tfas_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tfas_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  genes = Tfas_genes_multi[Tfas_genes_multi$chrom == chrom,]
  for (n in 1:n){
    seq = seq(starts_win[n], starts_win[n]+999999)
    counts = 0
    for (i in 1:nrow(genes)){
      in_window <- genes[i,3] %in% seq
      if (in_window == "TRUE"){
        counts = counts + 1
      } else {
        counts = counts
      }
    }
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, as.numeric(counts))
    Tfas_multi_gene_counts_per_mb_windows <- rbind(Tfas_multi_gene_counts_per_mb_windows, window_entry)
  }
}
#write.table(x = Tfas_multi_gene_counts_per_mb_windows, file = "Tfas_multicopy_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")
Tfas_multi_gene_counts_per_mb_windows <- read.table("Tfas_multicopy_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
summary(Tfas_multi_gene_counts_per_mb_windows)
colnames(Tfas_multi_gene_counts_per_mb_windows) <- c("chrom", "start_window", "end_window", "multicopy_gene_count")

## 4th track: REPETITIVE CONTENT

# Calculate repeat density
Tfas_te <- read.table("TE_table_fasciculata.for_circlize.txt")
colnames(Tfas_te) <- c("chrom", "class", "start", "end")
Tfas_te$length <- Tfas_te$end - Tfas_te$start

## Loop to create per 1 Mb window count of starting TEs
Tfas_te_counts_per_mb_windows <- data.frame()
for (i in 1:nrow(Tfas_chrom)){
  chrom = Tfas_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tfas_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  tes = Tfas_te[Tfas_te$chrom == chrom,]
  for (n in 1:n){
    seq = seq(starts_win[n], starts_win[n]+999999)
    counts = 0
    for (i in 1:nrow(tes)){
      in_window <- tes[i,3] %in% seq
      if (in_window == "TRUE"){
        counts = counts + 1
      } else {
        counts = counts
      }
    }
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, counts)
    Tfas_te_counts_per_mb_windows <- rbind(Tfas_te_counts_per_mb_windows, window_entry)
  }
}

## Loop to create per 1 Mb window the percentage of repetitive bases
Tfas_masked_percentage_per_mb_windows <- data.frame()
rest <- data.frame()
for (i in 1:nrow(Tfas_chrom)){
  chrom = Tfas_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tfas_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  chrom_tes = Tfas_te[Tfas_te$chrom == chrom,]
  for (n in 1:n){
    start <- starts_win[n]
    print(start)
    end <- starts_win[n]+1000000
    window_tes <- chrom_tes[chrom_tes$start > start & chrom_tes$start < end,]
    window_tes$length <- window_tes$end - window_tes$start
    if (dim(window_tes)[1] != 0){
      for (m in 1:nrow(window_tes)){
        if (m == 1){
        length = window_tes[m,5]
        } else {
          if (window_tes[m,3] >= window_tes[m-1,3] & window_tes[m,3] <= window_tes[m-1,4]
            & window_tes[m,4] > window_tes[m-1,4]){
            overlap <- window_tes[m-1,4] - window_tes[m,3]
            total <- length + window_tes[m,5]
            length = total - overlap
            } else if (window_tes[m,3] >= window_tes[m-1,3] & window_tes[m,4] <= window_tes[m-1,4]){
              length = length
            } else {
              length = length + window_tes[m,5]
        }
      }
      }
      rest_length <- rest$end - start
      if (dim(rest)[1] != 0){
      print("there's a rest!")
      te_length <- length + rest_length
    } else {
      te_length = length
    }
    percentage_te <- (te_length/1000000)*100
    } else {
      percentage_te <- 0
    }
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, as.numeric(percentage_te))
    Tfas_masked_percentage_per_mb_windows <- rbind(Tfas_masked_percentage_per_mb_windows, window_entry)
    rest <- window_tes[window_tes$end > end,]
  }
}
colnames(Tfas_masked_percentage_per_mb_windows) <- c("chrom", "start_window", "end_window", "percentage_te")
#write.table(x = Tfas_masked_percentage_per_mb_windows, file = "Tfas_TE_density_1Mb_windows.25scaffolds.txt")
Tfas_masked_percentage_per_mb_windows <- read.table("Tfas_TE_density_1Mb_windows.25scaffolds.txt")
Tfas_masked_percentage_per_mb_windows$percentage_te <- as.numeric(Tfas_masked_percentage_per_mb_windows$percentage_te)


# Make plot
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

### Introducing te proportions calculated with python (above loops in R are wrong)
Tfas_masked1 <- read.table("TE_content_Tfas_per1MBwindow_py.txt", header = T)
Tfas_masked_percentage_per_mb_windows <- Tfas_masked1
Tfas_masked2 <- read.table("Tfas_TE_density_1Mb_windows.25scaffolds.txt", header = T)
Tfas_masked_percentage_per_mb_windows <- as.data.frame(cbind(Tfas_masked2$chrom, Tfas_masked1$start_window,
                                               Tfas_masked1$end_window, Tfas_masked1$perc_te))
colnames(Tfas_masked_percentage_per_mb_windows) <- c("chrom", "start_window", "end_window", "percentage_te")
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


### For T. leiboldiana

# Select T.lei genes and remove all species specific orthogroups
Tlei_genes <- subset(genes_table, grepl("^Tlei_v1",V1))
Tlei_genes <- Tlei_genes[!(Tlei_genes$V8 == 0 & Tlei_genes$V9 == 0),]

# Calculate gene density
Tlei_chrom <- chrom[chrom$species=="Tlei",]

Tlei_genes <- Tlei_genes[,c(1,2,3,4)]
colnames(Tlei_genes) <- c("gene_id", "chrom", "start", "end")

## Loop to create per 1 Mb window count of starting genes
Tlei_all_gene_counts_per_mb_windows <- data.frame()
for (i in 1:nrow(Tlei_chrom)){
  chrom = Tlei_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tlei_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  genes = Tlei_genes[Tlei_genes$chrom == chrom,]
  for (n in 1:n){
    seq = seq(starts_win[n], starts_win[n]+999999)
    counts = 0
    for (i in 1:nrow(genes)){
      in_window <- genes[i,3] %in% seq
      if (in_window == "TRUE"){
        counts = counts + 1
      } else {
        counts = counts
      }
    }
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, counts)
    Tlei_all_gene_counts_per_mb_windows <- rbind(Tlei_all_gene_counts_per_mb_windows, window_entry)
  }
}
write.table(x = Tlei_all_gene_counts_per_mb_windows, file = "Tlei_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")
Tlei_all_gene_counts_per_mb_windows <- read.table("Tlei_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
summary(Tlei_all_gene_counts_per_mb_windows)

## Loop to create proportion of genic bases per 1 Mb window
Tlei_all_gene_content_per_1_mb_windows <- data.frame()
rest <- data.frame()
for (i in 1:nrow(Tlei_chrom)){
  chrom = Tlei_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tlei_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  chrom_genes = Tlei_genes[Tlei_genes$chrom == chrom,]
  for (n in 1:n){
    start <- starts_win[n]
    end <- starts_win[n]+1000000
    window_genes <- chrom_genes[chrom_genes$start > start & chrom_genes$start < end,]
    window_genes$length <- window_genes$end - window_genes$start
    genic_length <- sum(window_genes$length)
    rest_length <- rest$end - start
    if (dim(rest)[1] != 0){
      genic_length <- genic_length + rest_length
    }
    percentage_genic <- (genic_length/1000000)*100
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, as.numeric(percentage_genic))
    Tlei_all_gene_content_per_1_mb_windows <- rbind(Tlei_all_gene_content_per_1_mb_windows, window_entry)
    rest <- window_genes[window_genes$end > end,]
  }
}
write.table(Tlei_all_gene_content_per_1_mb_windows, file = "Tlei_gene_density_1Mb_windows.all_orthogroups_25_scaffolds.txt")
Tlei_all_gene_content_per_1_mb_windows <- read.table("Tlei_gene_density_1Mb_windows.all_orthogroups_25_scaffolds.txt")

## 2nd track: GENE DENSITY - 1:1 ORTHOGROUPS

# Select T.fas genes with 1:1 relationship
Tlei_genes_one <- Tlei_genes[(Tlei_genes$V9 == 1 & Tlei_genes$V10 == 1),]

Tlei_chrom <- chrom[chrom$species=="Tlei",]
Tlei_genes_one <- Tlei_genes_one[,c(1,2,3,4)]
colnames(Tlei_genes_one) <- c("gene_id", "chrom", "start", "end")

## Loop to create per 1 Mb window count of starting genes
Tlei_one_gene_counts_per_mb_windows <- data.frame()
for (i in 1:nrow(Tlei_chrom)){
  chrom = Tlei_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tlei_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  genes = Tlei_genes_one[Tlei_genes_one$chrom == chrom,]
  for (n in 1:n){
    start <- starts_win[n]
    end <- starts_win[n]+1000000
    counts <- nrow(genes[genes$start > start & genes$start < end,])
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, counts)
    Tlei_one_gene_counts_per_mb_windows <- rbind(Tlei_one_gene_counts_per_mb_windows, window_entry)
  }
}
write.table(x = Tlei_one_gene_counts_per_mb_windows, file = "Tlei_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")
Tlei_one_gene_counts_per_mb_windows <- read.table("Tlei_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
colnames(Tlei_one_gene_counts_per_mb_windows) <- c("chrom", "start_window", "end_window", "gene_count")
summary(Tlei_one_gene_counts_per_mb_windows)

## Loop to create proportion of genic bases per 1 Mb window
Tlei_one_gene_content_per_1_mb_windows <- data.frame()
rest <- data.frame()
for (i in 1:nrow(Tlei_chrom)){
  chrom = Tlei_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tlei_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  chrom_genes = Tlei_genes_one[Tlei_genes_one$chrom == chrom,]
  for (n in 1:n){
    start <- starts_win[n]
    end <- starts_win[n]+1000000
    window_genes <- chrom_genes[chrom_genes$start > start & chrom_genes$start < end,]
    window_genes$length <- window_genes$end - window_genes$start
    genic_length <- sum(window_genes$length)
    print(genic_length)
    rest_length <- rest$end - start
    if (dim(rest)[1] != 0){
      print("There's a rest! New length:")
      genic_length <- genic_length + rest_length
      print(genic_length)
    }
    percentage_genic <- (genic_length/1000000)*100
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, percentage_genic)
    Tlei_one_gene_content_per_1_mb_windows <- rbind(Tlei_one_gene_content_per_1_mb_windows, window_entry)
    rest <- window_genes[window_genes$end > end,]
  }
}
write.table(Tlei_one_gene_content_per_1_mb_windows, file = "Tlei_gene_density_1Mb_windows.one-to-one_orthogroups_25_scaffolds.txt")
Tlei_one_gene_content_per_1_mb_windows <- read.table("Tlei_gene_density_1Mb_windows.one-to-one_orthogroups_25_scaffolds.txt", header = T)
summary(Tlei_one_gene_content_per_1_mb_windows)

## 3nd track: All multi-copy genes in Tfas

Tlei_genes_multi <- Tlei_genes[!(Tlei_genes$V9 == 1),]
Tlei_genes_multi <- Tlei_genes_multi[,c(1,2,3,4)]
colnames(Tlei_genes_multi) <- c("gene_id", "chrom", "start", "end")

## Loop to create per 1 Mb window count of starting genes
Tlei_multi_gene_counts_per_mb_windows <- data.frame()
for (i in 1:nrow(Tlei_chrom)){
  chrom = Tlei_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tlei_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  genes = Tlei_genes_multi[Tlei_genes_multi$chrom == chrom,]
  for (n in 1:n){
    start <- starts_win[n]
    end <- starts_win[n]+1000000
    counts <- nrow(genes[genes$start > start & genes$start < end,])
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, counts)
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, as.numeric(counts))
    Tlei_multi_gene_counts_per_mb_windows <- rbind(Tlei_multi_gene_counts_per_mb_windows, window_entry)
  }
}
write.table(x = Tlei_multi_gene_counts_per_mb_windows, file = "Tlei_multicopy_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")
Tlei_multi_gene_counts_per_mb_windows <- read.table("Tlei_multicopy_counts_1Mb_windows_curated_ogs.25scaffolds.txt", header = T)
summary(Tlei_multi_gene_counts_per_mb_windows)
colnames(Tlei_multi_gene_counts_per_mb_windows) <- c("chrom", "start_window", "end_window", "multicopy_gene_count")

## 4th track: REPETITIVE CONTENT

Tlei_te <- read.table("TE_table_leiboldiana.for_circlize.txt")
colnames(Tlei_te) <- c("chrom", "class", "start", "end")
Tlei_te$length <- Tlei_te$end - Tlei_te$start

## Loop to create per 1 Mb window the percentage of repetitive bases
Tlei_masked_percentage_per_mb_windows <- data.frame()
rest <- data.frame()
for (i in 1:nrow(Tlei_chrom)){
  chrom = Tlei_chrom[i,1]
  print(chrom)
  starts_win <- seq(1, Tlei_chrom[i,2], by = 1000000)
  n <- length(starts_win)
  chrom_tes = Tlei_te[Tlei_te$chrom == chrom,]
  for (n in 1:n){
    start <- starts_win[n]
    print(start)
    end <- starts_win[n]+1000000
    window_tes <- chrom_tes[chrom_tes$start > start & chrom_tes$start < end,]
    window_tes$length <- window_tes$end - window_tes$start
    if (dim(window_tes)[1] != 0){
      for (m in 1:nrow(window_tes)){
        if (m == 1){
          length = window_tes[m,5]
        } else {
          if (window_tes[m,3] >= window_tes[m-1,3] & window_tes[m,3] <= window_tes[m-1,4]
              & window_tes[m,4] > window_tes[m-1,4]){
            overlap <- window_tes[m-1,4] - window_tes[m,3]
            total <- length + window_tes[m,5]
            length = total - overlap
          } else if (window_tes[m,3] >= window_tes[m-1,3] & window_tes[m,4] <= window_tes[m-1,4]){
            length = length
          } else {
            length = length + window_tes[m,5]
          }
        }
      }
      rest_length <- rest$end - start
      if (dim(rest)[1] != 0){
        print("there's a rest!")
        te_length <- length + rest_length
      } else {
        te_length = length
      }
      percentage_te <- (te_length/1000000)*100
    } else {
      percentage_te <- 0
    }
    window_entry <- cbind(chrom, starts_win[n], starts_win[n]+999999, as.numeric(percentage_te))
    Tlei_masked_percentage_per_mb_windows <- rbind(Tlei_masked_percentage_per_mb_windows, window_entry)
    rest <- window_tes[window_tes$end > end,]
  }
}
colnames(Tlei_masked_percentage_per_mb_windows) <- c("chrom", "start_window", "end_window", "percentage_te")
write.table(x = Tlei_masked_percentage_per_mb_windows, file = "Tlei_TE_density_1Mb_windows.25scaffolds.txt")
Tlei_masked_percentage_per_mb_windows <- read.table("Tlei_TE_density_1Mb_windows.25scaffolds.txt")
Tlei_masked_percentage_per_mb_windows$percentage_te <- as.numeric(Tlei_masked_percentage_per_mb_windows$percentage_te)

# Make plot

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

Tlei_masked1 <- read.table("TE_content_Tlei_per1MBwindow_py.txt", header = T)
Tlei_masked2 <- read.table("Tlei_TE_density_1Mb_windows.25scaffolds.txt", header = T)
Tlei_masked_percentage_per_mb_windows <- as.data.frame(cbind(Tlei_masked2$chrom, Tlei_masked1$start_window,
                                                             Tlei_masked1$end_window, Tlei_masked1$perc_te))
colnames(Tlei_masked_percentage_per_mb_windows) <- c("chrom", "start_window", "end_window", "percentage_te")
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
