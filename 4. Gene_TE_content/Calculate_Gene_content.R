## Set working directory and load needed packages
setwd("")

## Load a list of chromosomes for both genomes
chrom <- read.table("chromosomes_Tfas_Tlei.txt")
colnames(chrom) <- c("name", "size")
chrom$species <- c(rep(c("Tfas"), 25), rep(c("Tlei"), 26))

## Load in list of curated orthologous genes
genes_table <- read.delim("../orthology_Tfas_Tlei/orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt", sep = "\t", header = F)

#-----------------------------T.FASCICULATA------------------------------#

#1st track: GENE DENSITY - ALL GENES

# Select T.fas genes and remove all species specific orthogroups
Tfas_genes <- subset(genes_table, grepl("^Tfasc.v1",V1))
Tfas_genes <- Tfas_genes[!(Tfas_genes$V8 == 0 & Tfas_genes$V10 == 0),]
# Calculate gene count per window
Tfas_chrom <- chrom[chrom$species=="Tfas",]
Tfas_genes <- Tfas_genes[,c(1,2,3,4)]
colnames(Tfas_genes) <- c("gene_id", "chrom", "start", "end")

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
write.table(x = counts_per_mb_windows, file = "Tfas_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")


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
write.table(Tfas_all_gene_content_per_1_mb_windows, file = "Tfas_gene_density_1Mb_windows.all_orthogroups_25_scaffolds.txt")


## 2nd track: GENE DENSITY - 1:1 ORTHOGROUPS

# Select T.fas genes with 1:1 relationship
Tfas_genes_one <- Tfas_genes[(Tfas_genes$V9 == 1 & Tfas_genes$V10 == 1),]

Tfas_chrom <- chrom[chrom$species=="Tfas",]
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
write.table(x = Tfas_one_gene_counts_per_mb_windows, file = "Tfas_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")


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
write.table(Tfas_one_gene_content_per_1_mb_windows, file = "Tfas_gene_density_1Mb_windows.one-to-one_orthogroups_25_scaffolds.txt")


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
colnames(Tfas_multi_gene_counts_per_mb_windows) <- c("chrom", "start_window", "end_window", "multicopy_gene_count")
write.table(x = Tfas_multi_gene_counts_per_mb_windows, file = "Tfas_multicopy_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")



#-------------------------------T.LEIBOLDIANA------------------------------#

#1st track: GENE DENSITY - ALL GENES

# Select T.lei genes and remove all species specific orthogroups
Tlei_genes <- subset(genes_table, grepl("^Tlei_v1",V1))
Tlei_genes <- Tlei_genes[!(Tlei_genes$V8 == 0 & Tlei_genes$V9 == 0),]

# Loop to create per 1 Mb window count of starting genes
Tlei_chrom <- chrom[chrom$species=="Tlei",]
Tlei_genes <- Tlei_genes[,c(1,2,3,4)]
colnames(Tlei_genes) <- c("gene_id", "chrom", "start", "end")
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


## 2nd track: GENE DENSITY - 1:1 ORTHOGROUPS

# Select T.lei genes with 1:1 relationship
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
colnames(Tlei_one_gene_counts_per_mb_windows) <- c("chrom", "start_window", "end_window", "gene_count")
write.table(x = Tlei_one_gene_counts_per_mb_windows, file = "Tlei_one-to-one_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")

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

## 3nd track: MULTI-COPY genes

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
colnames(Tlei_multi_gene_counts_per_mb_windows) <- c("chrom", "start_window", "end_window", "multicopy_gene_count")
write.table(x = Tlei_multi_gene_counts_per_mb_windows, file = "Tlei_multicopy_gene_counts_1Mb_windows_curated_ogs.25scaffolds.txt")
