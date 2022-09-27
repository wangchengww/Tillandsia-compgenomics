#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2", "reshape2", "stringr","grid", "gridExtra", "matrixStats", "cowplot", "forcats")

# Extract legend to plot legend for multiple plots
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Load data
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/III. Expression figure/B-Expression-curves-multicopy/')
counts_Tfas <- read.table("counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)
counts_Tlei <- read.table("counts.Tfas_Tlei_6_timepoints.exons.toTLEI.sum.normalized-cpm.EdgeR.txt", header = T, row.names = 1)

genes <- scan("OG0001504_PEPC_gene_list.txt", character(), quote = "")

module_counts_Tfas <- subset(counts_Tfas, rownames(counts_Tfas) %in% genes)
module_counts_Tfas <- module_counts_Tfas[, c(1:36)]
module_counts_Tlei <- subset(counts_Tlei, rownames(counts_Tlei) %in% genes)
module_counts_Tlei <- module_counts_Tlei[, c(37:72)]
mod_size_Tfas <- nrow(module_counts_Tfas)
mod_size_Tlei <- nrow(module_counts_Tlei)

# Melt the count data
module_counts_Tfas$gene_id <- row.names(module_counts_Tfas)
module_counts_Tfas_m <- melt(module_counts_Tfas, id.vars = "gene_id")
module_counts_Tfas_m$time <- c(rep("N+9", mod_size_Tfas), rep("D+1",mod_size_Tfas), rep("D+5",mod_size_Tfas),
                               rep("D+9",mod_size_Tfas), rep("N+1",mod_size_Tfas), rep("N+5",mod_size_Tfas))
module_counts_Tfas_m$sample <- c(rep("A", as.numeric(mod_size_Tfas)*6), rep("B",as.numeric(mod_size_Tfas)*6),
                                 rep("C",as.numeric(mod_size_Tfas)*6), rep("D",as.numeric(mod_size_Tfas)*6),
                                 rep("E",as.numeric(mod_size_Tfas)*6), rep("F",as.numeric(mod_size_Tfas)*6))
colnames(module_counts_Tfas_m) <- c("gene_id", "id", "count", "time", "sample")

module_counts_Tlei$gene_id <- row.names(module_counts_Tlei)
module_counts_Tlei_m <- melt(module_counts_Tlei, id.vars = "gene_id")
module_counts_Tlei_m$time <- c(rep("N+9", mod_size_Tlei), rep("D+1",mod_size_Tlei), rep("D+5",mod_size_Tlei),
                               rep("D+9",mod_size_Tlei), rep("N+1",mod_size_Tlei), rep("N+5",mod_size_Tlei))
module_counts_Tlei_m$sample <- c(rep("A", as.numeric(mod_size_Tlei)*6), rep("B",as.numeric(mod_size_Tlei)*6),
                                 rep("C",as.numeric(mod_size_Tlei)*6), rep("D",as.numeric(mod_size_Tlei)*6),
                                 rep("E",as.numeric(mod_size_Tlei)*6), rep("F",as.numeric(mod_size_Tlei)*6))
colnames(module_counts_Tlei_m) <- c("gene_id", "id", "count", "time", "sample")

# Calculate mean counts per time point
mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tfas) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tlei) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")

sd_count_Tfas <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(sd_count_Tfas) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
sd_count_Tfas <- as.data.frame(sd_count_Tfas)
sd_count_Tlei <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(sd_count_Tlei) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
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
total_Tfas$species <- "T.fasciculata"
total_Tlei <- merge(mean_count_Tlei_m, sd_count_Tlei_m, by=c("gene_id", "time"))
total_Tlei$species <- "T.leiboldiana"

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

# Make plot 

total <- rbind(total_Tfas, total_Tlei)
total$order <- rep(c(4,5,6,1,2,3), 4)
nb.cols <- 4
mycolors <- brewer.pal(6, "Paired")
p1 <- ggplot(total, aes(x=fct_reorder(time, order), y=mean_count, group = gene_id)) +
  geom_point(aes(color=species)) +
  geom_line(aes(color=species, linetype = species), lwd = 1, ) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  #geom_vline(xintercept = 1.75, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, color = species), width=.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values=c("#3fc1c0", "#1a759f")) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "1) OG0001504: Phosphenolcarboxylase (PEPC) (2:2)") + theme(
    plot.title = element_text(size = 11)) +
  theme_bw()

pdf("Expression_curve.OG0001504_PEPC.27092022.pdf", width = 9, height = 6)
p1
dev.off()
###################### NEXT PLOT #########################

genes <- scan("OG0000580_Pyrophosphate−−fructose-6−phosphate-1−phosphotransferase-PFP_genelist.txt", character(), quote = "")

module_counts_Tfas <- subset(counts_Tfas, rownames(counts_Tfas) %in% genes)
module_counts_Tfas <- module_counts_Tfas[, c(1:36)]
module_counts_Tlei <- subset(counts_Tlei, rownames(counts_Tlei) %in% genes)
module_counts_Tlei <- module_counts_Tlei[, c(37:72)]
mod_size_Tfas <- nrow(module_counts_Tfas)
mod_size_Tlei <- nrow(module_counts_Tlei)

# Melt the count data
module_counts_Tfas$gene_id <- row.names(module_counts_Tfas)
module_counts_Tfas_m <- melt(module_counts_Tfas, id.vars = "gene_id")
module_counts_Tfas_m$time <- c(rep("N+9", mod_size_Tfas), rep("D+1",mod_size_Tfas), rep("D+5",mod_size_Tfas),
                               rep("D+9",mod_size_Tfas), rep("N+1",mod_size_Tfas), rep("N+5",mod_size_Tfas))
module_counts_Tfas_m$sample <- c(rep("A", as.numeric(mod_size_Tfas)*6), rep("B",as.numeric(mod_size_Tfas)*6),
                                 rep("C",as.numeric(mod_size_Tfas)*6), rep("D",as.numeric(mod_size_Tfas)*6),
                                 rep("E",as.numeric(mod_size_Tfas)*6), rep("F",as.numeric(mod_size_Tfas)*6))
colnames(module_counts_Tfas_m) <- c("gene_id", "id", "count", "time", "sample")

module_counts_Tlei$gene_id <- row.names(module_counts_Tlei)
module_counts_Tlei_m <- melt(module_counts_Tlei, id.vars = "gene_id")
module_counts_Tlei_m$time <- c(rep("N+9", mod_size_Tlei), rep("D+1",mod_size_Tlei), rep("D+5",mod_size_Tlei),
                               rep("D+9",mod_size_Tlei), rep("N+1",mod_size_Tlei), rep("N+5",mod_size_Tlei))
module_counts_Tlei_m$sample <- c(rep("A", as.numeric(mod_size_Tlei)*6), rep("B",as.numeric(mod_size_Tlei)*6),
                                 rep("C",as.numeric(mod_size_Tlei)*6), rep("D",as.numeric(mod_size_Tlei)*6),
                                 rep("E",as.numeric(mod_size_Tlei)*6), rep("F",as.numeric(mod_size_Tlei)*6))
colnames(module_counts_Tlei_m) <- c("gene_id", "id", "count", "time", "sample")

# Calculate mean counts per time point
mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tfas) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tlei) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")

sd_count_Tfas <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(sd_count_Tfas) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
sd_count_Tfas <- as.data.frame(sd_count_Tfas)
sd_count_Tlei <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(sd_count_Tlei) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
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
total_Tfas$species <- "T.fasciculata"
total_Tlei <- merge(mean_count_Tlei_m, sd_count_Tlei_m, by=c("gene_id", "time"))
total_Tlei$species <- "T.leiboldiana"

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

# Make plot 
total <- rbind(total_Tfas, total_Tlei)
total$order <- rep(c(4,5,6,1,2,3), 5)
p2 <- ggplot(total, aes(x=fct_reorder(time, order), y=mean_count, group = gene_id)) +
  geom_point(aes(color=species)) +
  geom_line(aes(color=species, linetype = species), lwd = 1, ) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, color = species), width=.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values=c("#a4133c", "#f6aa1c")) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "2) OG0000580: Pyrophosphate--fructose 6-phosphate 1-phosphotransferase (PFP) (3:2)") + theme(
    plot.title = element_text(size = 11)) +
  theme_bw()

pdf("Expression_curve.OG0000580_PFP.27092022.pdf", width = 9, height = 6)
p2
dev.off()
###################### NEXT PLOT ####################
genes <- scan("OG0005044_XAP5-CIRCADIAN-TIMEKEEPER_genelist.txt", character(), quote = "")

module_counts_Tfas <- subset(counts_Tfas, rownames(counts_Tfas) %in% genes)
module_counts_Tfas <- module_counts_Tfas[, c(1:36)]
module_counts_Tlei <- subset(counts_Tlei, rownames(counts_Tlei) %in% genes)
module_counts_Tlei <- module_counts_Tlei[, c(37:72)]
mod_size_Tfas <- nrow(module_counts_Tfas)
mod_size_Tlei <- nrow(module_counts_Tlei)

# Melt the count data
module_counts_Tfas$gene_id <- row.names(module_counts_Tfas)
module_counts_Tfas_m <- melt(module_counts_Tfas, id.vars = "gene_id")
module_counts_Tfas_m$time <- c(rep("N+9", mod_size_Tfas), rep("D+1",mod_size_Tfas), rep("D+5",mod_size_Tfas),
                               rep("D+9",mod_size_Tfas), rep("N+1",mod_size_Tfas), rep("N+5",mod_size_Tfas))
module_counts_Tfas_m$sample <- c(rep("A", as.numeric(mod_size_Tfas)*6), rep("B",as.numeric(mod_size_Tfas)*6),
                                 rep("C",as.numeric(mod_size_Tfas)*6), rep("D",as.numeric(mod_size_Tfas)*6),
                                 rep("E",as.numeric(mod_size_Tfas)*6), rep("F",as.numeric(mod_size_Tfas)*6))
colnames(module_counts_Tfas_m) <- c("gene_id", "id", "count", "time", "sample")

module_counts_Tlei$gene_id <- row.names(module_counts_Tlei)
module_counts_Tlei_m <- melt(module_counts_Tlei, id.vars = "gene_id")
module_counts_Tlei_m$time <- c(rep("N+9", mod_size_Tlei), rep("D+1",mod_size_Tlei), rep("D+5",mod_size_Tlei),
                               rep("D+9",mod_size_Tlei), rep("N+1",mod_size_Tlei), rep("N+5",mod_size_Tlei))
module_counts_Tlei_m$sample <- c(rep("A", as.numeric(mod_size_Tlei)*6), rep("B",as.numeric(mod_size_Tlei)*6),
                                 rep("C",as.numeric(mod_size_Tlei)*6), rep("D",as.numeric(mod_size_Tlei)*6),
                                 rep("E",as.numeric(mod_size_Tlei)*6), rep("F",as.numeric(mod_size_Tlei)*6))
colnames(module_counts_Tlei_m) <- c("gene_id", "id", "count", "time", "sample")

# Calculate mean counts per time point
mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tfas) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
mean_count_Tlei <- t(mean_count_Tlei)
colnames(mean_count_Tlei) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
rownames(mean_count_Tlei) <- rownames(module_counts_Tlei)
mean_count_Tlei <- as.data.frame(mean_count_Tlei)

sd_count_Tfas <- sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(sd_count_Tfas) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
sd_count_Tfas <- as.data.frame(sd_count_Tfas)
sd_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowSds(as.matrix(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)]))))
sd_count_Tlei <- t(sd_count_Tlei)
colnames(sd_count_Tlei) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
rownames(sd_count_Tlei) <- rownames(module_counts_Tlei)
sd_count_Tlei <- as.data.frame(sd_count_Tlei)

mean_count_Tfas$gene_id <- row.names(mean_count_Tfas)
mean_count_Tlei$gene_id <- row.names(mean_count_Tlei)
sd_count_Tfas$gene_id <- row.names(mean_count_Tfas)
sd_count_Tlei$gene_id <- row.names(mean_count_Tlei)

mean_count_Tfas_m <- melt(mean_count_Tfas, id.vars = "gene_id")
colnames(mean_count_Tfas_m) <- c("gene_id", "time", "mean_count")
mean_count_Tfas_m$time <- as.character(mean_count_Tfas_m$time)
mean_count_Tlei_m <- melt(mean_count_Tlei, id.vars = "gene_id")
colnames(mean_count_Tlei_m) <- c("gene_id", "time", "mean_count")
mean_count_Tlei_m$time <- as.character(mean_count_Tlei_m$time)

sd_count_Tfas_m <- melt(sd_count_Tfas, id.vars = "gene_id")
colnames(sd_count_Tfas_m) <- c("gene_id", "time", "sd")
sd_count_Tfas_m$time <- as.character(sd_count_Tfas_m$time)
sd_count_Tlei_m <- melt(sd_count_Tlei, id.vars = "gene_id")
colnames(sd_count_Tlei_m) <- c("gene_id", "time", "sd")
sd_count_Tlei_m$time <- as.character(sd_count_Tlei_m$time)


total_Tfas <- merge(mean_count_Tfas_m, sd_count_Tfas_m, by=c("gene_id", "time"))
total_Tfas$species <- "T.fasciculata"
total_Tlei <- merge(mean_count_Tlei_m, sd_count_Tlei_m, by=c("gene_id", "time"))
total_Tlei$species <- "T.leiboldiana"

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

# Make plot
total <- rbind(total_Tfas, total_Tlei)
total$order <- rep(c(4,5,6,1,2,3), 3)
p3 <- ggplot(total, aes(x=fct_reorder(time, order), y=mean_count, group = gene_id)) +
  geom_point(aes(color=species)) +
  geom_line(aes(color=species, linetype = species), lwd = 1, ) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  geom_errorbar(aes(ymin=mean_count-sd, ymax=mean_count+sd, color = species), width=.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values=c("#c77dff", "#7b2cbf")) +
  ylim(c(ymin, ymax)) +
  ylab("Counts Per Million (CPM)") +
  xlab("Time") +
  ggtitle(label = "3) OG0005044: XAP5 CIRCADIAN TIMEKEEPER (2:1)") + theme(
    plot.title = element_text(size = 11)) +
  theme_bw()

pdf("Expression_curve.OG0005044_XAP5.27092022.pdf", width = 9, height = 6)
p3
dev.off()

pdf("Expression_curve_MULTICOPY_withTlei.pdf", width = 8, height = 12)
grid.arrange(p1, p2, p3, ncol = 1, nrow = 3)
dev.off()
