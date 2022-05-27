#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2", "reshape2", "stringr","grid", "gridExtra", "matrixStats", "cowplot")

# Extract legend to plot legend for multiple plots
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Load data
setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/B-Expression-curves-multicopy/')
counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)
genes <- scan("OG0001504_PEPC_gene_list.txt", character(), quote = "")

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
  scale_color_manual(values=c("#3fc1c0", "#1a759f")) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "T. fasciculata") + theme(
    plot.title = element_text(size = 11)) +
  theme_bw()

# In Tlei
counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.toTLEI.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)
genes <- scan("OG0001504_PEPC_gene_list.txt", character(), quote = "")

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
  scale_color_manual(values=c("#3fc1c0", "#1a759f")) +
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
plot1 <- grid.arrange(p1, p2, nrow = 1, widths=c(2.3, 2.3, 1), 
                      top = grid::textGrob("1) OG0001504: Phosphenolcarboxylase (PEPC)", x = 0.035, hjust = 0, gp=gpar(fontsize=13)))

###################### NEXT PLOT #########################

genes <- scan("OG0000580_PFP_gene_list_withTlei.txt", character(), quote = "")
counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)

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
p3 <- ggplot(total, aes(x=time, y=mean_count, group = interaction(gene_id,species))) +
  geom_point(aes(color=gene_id)) +
  geom_line(aes(color=gene_id, linetype = species), lwd = 1) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, color = gene_id), width=.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values=c('#da9e3f','#a5323d', '#327c5e')) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "T. fasciculata") + theme(plot.title = element_text(size = 11))+
  theme_bw()

# In Tlei
counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.toTLEI.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)

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
min_Tlei <- min(total_Tlei$mean_count+total_Tlei$sd)

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
p4 <- ggplot(total, aes(x=time, y=mean_count, group = interaction(gene_id,species))) +
  geom_point(aes(color=gene_id)) +
  geom_line(aes(color=gene_id, linetype = species), lwd = 1) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, color = gene_id), width=.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values=c('#da9e3f','#a5323d', '#327c5e')) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "T. leiboldiana") + theme(plot.title = element_text(size = 11))+
  theme_bw()

# Extract legend
#legend1 <- get_legend(p1)
#Remove the legend from the box plot
#p1 <- p1 + theme(legend.position="none")
#p2 <- p2 + theme(legend.position="none")
# Arrange ggplot2 graphs with a specific width
plot2 <- grid.arrange(p3, p4, nrow = 1, widths=c(2.3, 2.3, 1), 
                      top = grid::textGrob("2) OG0000580: Pyrophosphate--fructose 6-phosphate 1-phosphotransferase (PFP)", x = 0.035, hjust = 0, gp=gpar(fontsize=13)))

###################### NEXT PLOT ####################
genes <- scan("OG0005044_XCT_list_genes_withTlei.txt", character(), quote = "")

counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)

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
p5 <- ggplot(total, aes(x=time, y=mean_count, group = interaction(gene_id,species))) +
  geom_point(aes(color=gene_id)) +
  geom_line(aes(color=gene_id, linetype = species), lwd = 1) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, color = gene_id), width=.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values=c("#c77dff", "#7b2cbf")) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "T. fasciculata") + theme(plot.title = element_text(size = 11))+
  theme_bw()

# In Tlei
counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.toTLEI.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)

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
p6 <- ggplot(total, aes(x=time, y=mean_count, group = interaction(gene_id,species))) +
  geom_point(aes(color=gene_id)) +
  geom_line(aes(color=gene_id, linetype = species), lwd = 1) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, color = gene_id), width=.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values=c("#c77dff", "#7b2cbf")) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "T. leiboldiana") + theme(plot.title = element_text(size = 11))+
  theme_bw()

# Arrange ggplot2 graphs with a specific width
plot3 <- grid.arrange(p5, p6, nrow = 1, widths=c(2.3, 2.3, 1), 
                      top = grid::textGrob("3) OG0005044: XAP5 CIRCADIAN TIMEKEEPER", x = 0.035, hjust = 0, gp=gpar(fontsize=13)))

grid.arrange(plot1, plot2, plot3, ncol = 1, nrow = 3)

pdf("Expression_curve_MULTICOPY.pdf", width = 10, height = 16)
grid.arrange(plot1, plot2, plot3, ncol = 1, nrow = 3)
dev.off()
