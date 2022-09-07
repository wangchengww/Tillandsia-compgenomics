#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2", "reshape2", "stringr","grid", "gridExtra", 
               "matrixStats", "cowplot", "RColorBrewer")

# Extract legend to plot legend for multiple plots
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Obtain arguments from command line. 1 = counts in T. fasciculata, 2 = counts in
# T. leiboldiana, 3 = gene list for orthogroup, 4 = title of the plot, 5 = color palette
args <- commandArgs(trailingOnly = TRUE)

nb.cols <- 8
#mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb.cols)
mycolors <- brewer.pal(8, args[[4]])

# Load data
#setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/B-Expression-curves-multicopy/')
#counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)
counts <- read.table(args[[1]], header = T, row.names = 1)
#genes <- scan("OG0001504_PEPC_gene_list_withTlei.txt", character(), quote = "")
genes <- scan(args[[3]], character(), quote = "")
og <- str_split(args[[3]], "\\_|\\.")[[1]][1]
func <- str_split(args[[3]], "\\_|\\.")[[1]][2]

module_counts <- subset(counts, rownames(counts) %in% genes)
module_counts_Tfas <- module_counts[, c(1:36)]
module_counts_Tlei <- module_counts[, c(37:72)]
mod_size <- as.character(nrow(module_counts))

# Melt the count data
module_counts_Tfas$gene_id <- row.names(module_counts_Tfas)
module_counts_Tfas_m <- melt(module_counts_Tfas, id.vars = "gene_id")
module_counts_Tfas_m$time <- c(rep("0100", mod_size), rep("0500",mod_size), rep("0900",mod_size),
                               rep("1300",mod_size), rep("1700",mod_size), rep("2100",mod_size))
module_counts_Tfas_m$sample <- c(rep("A", as.numeric(mod_size)*6), rep("B",as.numeric(mod_size)*6),
                                 rep("C",as.numeric(mod_size)*6), rep("D",as.numeric(mod_size)*6),
                                 rep("E",as.numeric(mod_size)*6), rep("F",as.numeric(mod_size)*6))
colnames(module_counts_Tfas_m) <- c("gene_id", "id", "count", "time", "sample")

module_counts_Tlei$gene_id <- row.names(module_counts_Tlei)
module_counts_Tlei_m <- melt(module_counts_Tlei, id.vars = "gene_id")
module_counts_Tlei_m$time <- c(rep("0100", mod_size), rep("0500",mod_size), rep("0900",mod_size),
                               rep("1300",mod_size), rep("1700",mod_size), rep("2100",mod_size))
module_counts_Tlei_m$sample <- c(rep("A", as.numeric(mod_size)*6), rep("B",as.numeric(mod_size)*6),
                                 rep("C",as.numeric(mod_size)*6), rep("D",as.numeric(mod_size)*6),
                                 rep("E",as.numeric(mod_size)*6), rep("F",as.numeric(mod_size)*6))
colnames(module_counts_Tlei_m) <- c("gene_id", "id", "count", "time", "sample")

# Calculate mean counts per time point
if (nrow(module_counts_Tfas) == 1){
  mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  mean_count_Tfas <- t(mean_count_Tfas)
  colnames(mean_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
  rownames(mean_count_Tfas) <- rownames(module_counts_Tfas)
  mean_count_Tfas <- as.data.frame(mean_count_Tfas)
  mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  mean_count_Tlei <- t(mean_count_Tlei)
  colnames(mean_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")
  rownames(mean_count_Tlei) <- rownames(module_counts_Tlei)
  mean_count_Tlei <- as.data.frame(mean_count_Tlei)
  
  sd_count_Tfas <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  sd_count_Tfas <- t(sd_count_Tfas)
  colnames(sd_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
  rownames(sd_count_Tfas) <- rownames(module_counts_Tfas)
  sd_count_Tfas <- as.data.frame(sd_count_Tfas)
  sd_count_Tlei <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  sd_count_Tlei <- t(sd_count_Tlei)
  colnames(sd_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")
  rownames(sd_count_Tlei) <- rownames(module_counts_Tlei)
  sd_count_Tlei <- as.data.frame(sd_count_Tlei)
} else {
  mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  colnames(mean_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
  mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  colnames(mean_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")
  
  sd_count_Tfas <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  colnames(sd_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
  sd_count_Tfas <- as.data.frame(sd_count_Tfas)
  sd_count_Tlei <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  colnames(sd_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")
  sd_count_Tlei <- as.data.frame(sd_count_Tlei)
}

mean_count_Tfas$gene_id <- row.names(mean_count_Tfas)
mean_count_Tlei$gene_id <- row.names(mean_count_Tlei)
sd_count_Tfas$gene_id <- row.names(mean_count_Tfas)
sd_count_Tlei$gene_id <- row.names(mean_count_Tlei)

mean_count_Tfas_m <- melt(mean_count_Tfas, id.vars = "gene_id")
colnames(mean_count_Tfas_m) <- c("gene_id", "time", "mean_count")
mean_count_Tlei_m <- melt(mean_count_Tlei, id.vars = "gene_id")
colnames(mean_count_Tlei_m) <- c("gene_id", "time", "mean_count")

sd_count_Tfas_m <- melt(sd_count_Tfas, id.vars = "gene_id")
colnames(sd_count_Tfas_m) <- c("gene_id", "time", "sd")
sd_count_Tlei_m <- melt(sd_count_Tlei, id.vars = "gene_id")
colnames(sd_count_Tlei_m) <- c("gene_id", "time", "sd")

total_Tfas <- merge(mean_count_Tfas_m, sd_count_Tfas_m, by=c("gene_id", "time"))
total_Tlei <- merge(mean_count_Tlei_m, sd_count_Tlei_m, by=c("gene_id", "time"))


# Calculate scale
max_Tfas <- max(total_Tfas$mean_count+total_Tfas$sd)
min_Tfas <- min(total_Tfas$mean_count-total_Tfas$sd)
max_Tlei <- max(total_Tlei$mean_count+total_Tlei$sd)
min_Tlei <- min(total_Tlei$mean_count-total_Tlei$sd)

if (max_Tfas > max_Tlei){
  ymax = max_Tfas + 1
} else {
  ymax = max_Tlei + 1
}

if (min_Tfas < min_Tlei){
  ymin = min_Tfas - 1
} else {
  ymin = min_Tlei - 1
}

# Highlight genes of interest
total_Tfas$species <- "T.fasciculata"
total_Tlei$species <- "T.leiboldiana"
total <- rbind(total_Tfas, total_Tlei)
p1 <- ggplot(total, aes(x=time, y=mean_count, group = interaction(gene_id,species))) +
  geom_point(aes(color=gene_id)) +
  geom_line(aes(color=gene_id, linetype = species), lwd = 1) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, color = gene_id), width=.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values=mycolors[3:(as.numeric(mod_size)+3)]) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "T. fasciculata") + theme(
    plot.title = element_text(size = 11)) +
  theme_bw()

# In Tlei
start_colours <- as.numeric(mod_size)+3
#counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.toTLEI.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)
counts <- read.table(args[[2]], header = T, row.names = 1)

module_counts <- subset(counts, rownames(counts) %in% genes)
module_counts_Tfas <- module_counts[, c(1:36)]
module_counts_Tlei <- module_counts[, c(37:72)]
mod_size <- as.character(nrow(module_counts))

# Melt the count data
module_counts_Tfas$gene_id <- row.names(module_counts_Tfas)
module_counts_Tfas_m <- melt(module_counts_Tfas, id.vars = "gene_id")
module_counts_Tfas_m$time <- c(rep("0100", mod_size), rep("0500",mod_size), rep("0900",mod_size),
                               rep("1300",mod_size), rep("1700",mod_size), rep("2100",mod_size))
module_counts_Tfas_m$sample <- c(rep("A", as.numeric(mod_size)*6), rep("B",as.numeric(mod_size)*6),
                                 rep("C",as.numeric(mod_size)*6), rep("D",as.numeric(mod_size)*6),
                                 rep("E",as.numeric(mod_size)*6), rep("F",as.numeric(mod_size)*6))
colnames(module_counts_Tfas_m) <- c("gene_id", "id", "count", "time", "sample")

module_counts_Tlei$gene_id <- row.names(module_counts_Tlei)
module_counts_Tlei_m <- melt(module_counts_Tlei, id.vars = "gene_id")
module_counts_Tlei_m$time <- c(rep("0100", mod_size), rep("0500",mod_size), rep("0900",mod_size),
                               rep("1300",mod_size), rep("1700",mod_size), rep("2100",mod_size))
module_counts_Tlei_m$sample <- c(rep("A", as.numeric(mod_size)*6), rep("B",as.numeric(mod_size)*6),
                                 rep("C",as.numeric(mod_size)*6), rep("D",as.numeric(mod_size)*6),
                                 rep("E",as.numeric(mod_size)*6), rep("F",as.numeric(mod_size)*6))
colnames(module_counts_Tlei_m) <- c("gene_id", "id", "count", "time", "sample")

# Calculate mean counts per time point
if (nrow(module_counts_Tfas) == 1){
  mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  mean_count_Tfas <- t(mean_count_Tfas)
  colnames(mean_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
  rownames(mean_count_Tfas) <- rownames(module_counts_Tfas)
  mean_count_Tfas <- as.data.frame(mean_count_Tfas)
  mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  mean_count_Tlei <- t(mean_count_Tlei)
  colnames(mean_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")
  rownames(mean_count_Tlei) <- rownames(module_counts_Tlei)
  mean_count_Tlei <- as.data.frame(mean_count_Tlei)
  
  sd_count_Tfas <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  sd_count_Tfas <- t(sd_count_Tfas)
  colnames(sd_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
  rownames(sd_count_Tfas) <- rownames(module_counts_Tfas)
  sd_count_Tfas <- as.data.frame(sd_count_Tfas)
  sd_count_Tlei <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  sd_count_Tlei <- t(sd_count_Tlei)
  colnames(sd_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")
  rownames(sd_count_Tlei) <- rownames(module_counts_Tlei)
  sd_count_Tlei <- as.data.frame(sd_count_Tlei)
} else {
  mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  colnames(mean_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
  mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  colnames(mean_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")
  
  sd_count_Tfas <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  colnames(sd_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
  sd_count_Tfas <- as.data.frame(sd_count_Tfas)
  sd_count_Tlei <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
  colnames(sd_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")
  sd_count_Tlei <- as.data.frame(sd_count_Tlei)
}

mean_count_Tfas$gene_id <- row.names(mean_count_Tfas)
mean_count_Tlei$gene_id <- row.names(mean_count_Tlei)
sd_count_Tfas$gene_id <- row.names(mean_count_Tfas)
sd_count_Tlei$gene_id <- row.names(mean_count_Tlei)

mean_count_Tfas_m <- melt(mean_count_Tfas, id.vars = "gene_id")
colnames(mean_count_Tfas_m) <- c("gene_id", "time", "mean_count")
mean_count_Tlei_m <- melt(mean_count_Tlei, id.vars = "gene_id")
colnames(mean_count_Tlei_m) <- c("gene_id", "time", "mean_count")

sd_count_Tfas_m <- melt(sd_count_Tfas, id.vars = "gene_id")
colnames(sd_count_Tfas_m) <- c("gene_id", "time", "sd")
sd_count_Tlei_m <- melt(sd_count_Tlei, id.vars = "gene_id")
colnames(sd_count_Tlei_m) <- c("gene_id", "time", "sd")

total_Tfas <- merge(mean_count_Tfas_m, sd_count_Tfas_m, by=c("gene_id", "time"))
total_Tlei <- merge(mean_count_Tlei_m, sd_count_Tlei_m, by=c("gene_id", "time"))

# Calculate scale
max_Tfas <- max(total_Tfas$mean_count+total_Tfas$sd)
min_Tfas <- min(total_Tfas$mean_count-total_Tfas$sd)
max_Tlei <- max(total_Tlei$mean_count+total_Tlei$sd)
min_Tlei <- min(total_Tlei$mean_count-total_Tlei$sd)

if (max_Tfas > max_Tlei){
  ymax = max_Tfas + 1
} else {
  ymax = max_Tlei + 1
}

if (min_Tfas < min_Tlei){
  ymin = min_Tfas - 1
} else {
  ymin = min_Tlei - 1
}

# Highlight genes of interest
total_Tfas$species <- "T.fasciculata"
total_Tlei$species <- "T.leiboldiana"
total <- rbind(total_Tfas, total_Tlei)
p2 <- ggplot(total, aes(x=time, y=mean_count, group = interaction(gene_id,species))) +
  geom_point(aes(color=gene_id)) +
  geom_line(aes(color=gene_id, linetype = species), lwd = 1) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, color = gene_id), width=.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values=mycolors[start_colours:(as.numeric(mod_size)+start_colours)]) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "T. leiboldiana") + theme(
    plot.title = element_text(size = 11)) +
  theme_bw()

# Extract legend
#legend1 <- get_legend(p1)
#Remove the legend from the box plot
#p1 <- p1 + theme(legend.position="none")
#p2 <- p2 + theme(legend.position="none")
# Arrange ggplot2 graphs with a specific width
pdf(paste0("Expression_curve_", og, "_", func, ".pdf"), width = 12, height = 6)
plot1 <- grid.arrange(p1, p2, nrow = 1, 
                      top = grid::textGrob(paste0(og, ": ", func), x = 0.035, hjust = 0, gp=gpar(fontsize=13)))
#plot1 <- grid.arrange(p1, p2, nrow = 1, widths=c(2.3, 2.3, 1), 
#                      top = grid::textGrob("wawa", x = 0.035, hjust = 0, gp=gpar(fontsize=13)))
