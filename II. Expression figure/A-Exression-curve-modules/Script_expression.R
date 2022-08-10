#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2","grid", "gridExtra")

setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/A-Exression-curve-modules')
tfas_cluster5 <- read.table('Mean_Expression_counts_T.fasciculata-cluster5.txt', header = T)
tlei_cluster5 <- read.table('Mean_Expression_counts_T.leiboldiana-cluster5.txt', header = T)

tfas_cluster1 <- read.table('Mean_Expression_counts_T.fasciculata-cluster1.txt', header = T)
tlei_cluster1 <- read.table('Mean_Expression_counts_T.leiboldiana-cluster1.txt', header = T)

tfas_cluster3 <- read.table('Mean_Expression_counts_T.fasciculata-cluster3.txt', header = T)
tlei_cluster3 <- read.table('Mean_Expression_counts_T.leiboldiana-cluster3.txt', header = T)


tfas_cluster5$time <- as.character(tfas_cluster5$time)
tfas_cluster5[tfas_cluster5$time == '100',]$time <- "0100"
tfas_cluster5[tfas_cluster5$time == '500',]$time <- "0500"
tfas_cluster5[tfas_cluster5$time == '900',]$time <- "0900"

tlei_cluster5$time <- as.character(tlei_cluster5$time)
tlei_cluster5[tlei_cluster5$time == '100',]$time <- "0100"
tlei_cluster5[tlei_cluster5$time == '500',]$time <- "0500"
tlei_cluster5[tlei_cluster5$time == '900',]$time <- "0900"

tfas_cluster1$time <- as.character(tfas_cluster1$time)
tfas_cluster1[tfas_cluster1$time == '100',]$time <- "0100"
tfas_cluster1[tfas_cluster1$time == '500',]$time <- "0500"
tfas_cluster1[tfas_cluster1$time == '900',]$time <- "0900"

tlei_cluster1$time <- as.character(tlei_cluster1$time)
tlei_cluster1[tlei_cluster1$time == '100',]$time <- "0100"
tlei_cluster1[tlei_cluster1$time == '500',]$time <- "0500"
tlei_cluster1[tlei_cluster1$time == '900',]$time <- "0900"

tfas_cluster3$time <- as.character(tfas_cluster3$time)
tfas_cluster3[tfas_cluster3$time == '100',]$time <- "0100"
tfas_cluster3[tfas_cluster3$time == '500',]$time <- "0500"
tfas_cluster3[tfas_cluster3$time == '900',]$time <- "0900"

tlei_cluster3$time <- as.character(tlei_cluster3$time)
tlei_cluster3[tlei_cluster3$time == '100',]$time <- "0100"
tlei_cluster3[tlei_cluster3$time == '500',]$time <- "0500"
tlei_cluster3[tlei_cluster3$time == '900',]$time <- "0900"

# Calculate scale
max_Tfas_cluster5 <- max(tfas_cluster5$count)
min_Tfas_cluster5 <- min(tfas_cluster5$count)
max_Tlei_cluster5 <- max(tlei_cluster5$count)
min_Tlei_cluster5 <- min(tlei_cluster5$count)

if (max_Tfas_cluster5 > max_Tlei_cluster5){
  ymax_cluster5 = max_Tfas_cluster5 + 1
} else {
  ymax_cluster5 = max_Tlei_cluster5 + 1
}

if (min_Tfas_cluster5 < min_Tlei_cluster5){
  ymin_cluster5 = min_Tfas_cluster5 - 1
} else {
  ymin_cluster5 = min_Tlei_cluster5 - 1
}

max_Tfas_cluster1 <- max(tfas_cluster1$count)
min_Tfas_cluster1 <- min(tfas_cluster1$count)
max_Tlei_cluster1 <- max(tlei_cluster1$count)
min_Tlei_cluster1 <- min(tlei_cluster1$count)

if (max_Tfas_cluster1 > max_Tlei_cluster1){
  ymax_cluster1 = max_Tfas_cluster1 + 1
} else {
  ymax_cluster1 = max_Tlei_cluster1 + 1
}

if (min_Tfas_cluster1 < min_Tlei_cluster1){
  ymin_cluster1 = min_Tfas_cluster1 - 1
} else {
  ymin_cluster1 = min_Tlei_cluster1 - 1
}

max_Tfas_cluster3 <- max(tfas_cluster3$count)
min_Tfas_cluster3 <- min(tfas_cluster3$count)
max_Tlei_cluster3 <- max(tlei_cluster3$count)
min_Tlei_cluster3 <- min(tlei_cluster3$count)

if (max_Tfas_cluster3 > max_Tlei_cluster3){
  ymax_cluster3 = max_Tfas_cluster3 + 1
} else {
  ymax_cluster3 = max_Tlei_cluster3 + 1
}

if (min_Tfas_cluster3 < min_Tlei_cluster3){
  ymin_cluster3 = min_Tfas_cluster3 - 1
} else {
  ymin_cluster3 = min_Tlei_cluster3 - 1
}
# Calculate total mean curve
medians_Tfas_cluster5 <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
                         count = c(median(tfas_cluster5[tfas_cluster5$time == '0100',]$count),
                                   median(tfas_cluster5[tfas_cluster5$time == '0500',]$count),
                                   median(tfas_cluster5[tfas_cluster5$time == '0900',]$count),
                                   median(tfas_cluster5[tfas_cluster5$time == '1300',]$count),
                                   median(tfas_cluster5[tfas_cluster5$time == '1700',]$count),
                                   median(tfas_cluster5[tfas_cluster5$time == '2100',]$count)))

medians_Tlei_cluster5 <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
                           count = c(median(tlei_cluster5[tlei_cluster5$time == '0100',]$count),
                                     median(tlei_cluster5[tlei_cluster5$time == '0500',]$count),
                                     median(tlei_cluster5[tlei_cluster5$time == '0900',]$count),
                                     median(tlei_cluster5[tlei_cluster5$time == '1300',]$count),
                                     median(tlei_cluster5[tlei_cluster5$time == '1700',]$count),
                                     median(tlei_cluster5[tlei_cluster5$time == '2100',]$count)))

medians_Tfas_cluster1 <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
                                    count = c(median(tfas_cluster1[tfas_cluster1$time == '0100',]$count),
                                              median(tfas_cluster1[tfas_cluster1$time == '0500',]$count),
                                              median(tfas_cluster1[tfas_cluster1$time == '0900',]$count),
                                              median(tfas_cluster1[tfas_cluster1$time == '1300',]$count),
                                              median(tfas_cluster1[tfas_cluster1$time == '1700',]$count),
                                              median(tfas_cluster1[tfas_cluster1$time == '2100',]$count)))

medians_Tlei_cluster1 <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
                                    count = c(median(tlei_cluster1[tlei_cluster1$time == '0100',]$count),
                                              median(tlei_cluster1[tlei_cluster1$time == '0500',]$count),
                                              median(tlei_cluster1[tlei_cluster1$time == '0900',]$count),
                                              median(tlei_cluster1[tlei_cluster1$time == '1300',]$count),
                                              median(tlei_cluster1[tlei_cluster1$time == '1700',]$count),
                                              median(tlei_cluster1[tlei_cluster1$time == '2100',]$count)))

medians_Tfas_cluster3 <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
                                    count = c(median(tfas_cluster3[tfas_cluster3$time == '0100',]$count),
                                              median(tfas_cluster3[tfas_cluster3$time == '0500',]$count),
                                              median(tfas_cluster3[tfas_cluster3$time == '0900',]$count),
                                              median(tfas_cluster3[tfas_cluster3$time == '1300',]$count),
                                              median(tfas_cluster3[tfas_cluster3$time == '1700',]$count),
                                              median(tfas_cluster3[tfas_cluster3$time == '2100',]$count)))

medians_Tlei_cluster3 <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
                                    count = c(median(tlei_cluster3[tlei_cluster3$time == '0100',]$count),
                                              median(tlei_cluster3[tlei_cluster3$time == '0500',]$count),
                                              median(tlei_cluster3[tlei_cluster3$time == '0900',]$count),
                                              median(tlei_cluster3[tlei_cluster3$time == '1300',]$count),
                                              median(tlei_cluster3[tlei_cluster3$time == '1700',]$count),
                                              median(tlei_cluster3[tlei_cluster3$time == '2100',]$count)))
# Highlight genes of interest
pdf("Expression_curve_PANEL_modules1_5_3.newcolors.pdf", width = 10, height = 16)
######################### CLUSTER 5 #############################
p1a <- ggplot(tfas_cluster5, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = tfas_cluster5[tfas_cluster5$gene_id == "Tfasc_v1.16311-RA",], aes(group = gene_id),
            size = .8, color = "#caa8f5") +
  geom_line(data = tfas_cluster5[tfas_cluster5$gene_id == "Tfasc_v1.03128-RA",], aes(group = gene_id),
            size = .8, color = "#83c5be") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(c(ymin_cluster5, ymax_cluster5)) +
  ggtitle("1A) cluster 5 - T. fasciculata") +
  ylab("Median-centered log(CPM)") +
  xlab("Time")

p1b <- ggplot(tfas_cluster5, aes(x=time, y=count, group = gene_id)) +
  geom_line(data = medians_Tfas_cluster5, aes(group = 1), size = 1, color = "black") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

p1 <- p1a + annotation_custom(ggplotGrob(p1b), xmin = 5, xmax = 6.5,
                              ymin = -9, ymax = -6)

p2a <- ggplot(tlei_cluster5, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = tlei_cluster5[tlei_cluster5$gene_id == "Tfasc_v1.16311-RA",], aes(group = gene_id),
            size = .8, color = "#9984d4") +
  geom_line(data = tlei_cluster5[tlei_cluster5$gene_id == "Tfasc_v1.03128-RA",], aes(group = gene_id),
            size = .8, color = "#006d77") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(ymin_cluster5, ymax_cluster5) +
  ggtitle("1B) cluster 5 - T. leiboldiana") +
  ylab("Median-centered log(CPM)") +
  xlab("Time")

p2b <- ggplot(tlei_cluster5, aes(x=time, y=count, group = gene_id)) +
  geom_line(data = medians_Tlei_cluster5, aes(group = 1), size = 1, color = "black") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  ylim(c(-1.35,-0.75))

p2 <- p2a + annotation_custom(ggplotGrob(p2b), xmin = 5, xmax = 6.5,
                              ymin = -9, ymax = -6)

######################### CLUSTER 1 #############################
p3a <- ggplot(tfas_cluster1, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = tfas_cluster1[tfas_cluster1$gene_id == "Tfasc_v1.09028-RA",], aes(group = gene_id),
            size = .8, color = "darkgreen") +
  geom_line(data = tfas_cluster1[tfas_cluster1$gene_id == "Tfasc_v1.19779-RA",], aes(group = gene_id),
            size = .8, color = "darkgreen") +
  geom_line(data = tfas_cluster1[tfas_cluster1$gene_id == "Tfasc_v1.08739-RA",], aes(group = gene_id),
            size = .8, color = "darkgreen") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(c(ymin_cluster1, ymax_cluster1)) +
  ggtitle("2A) cluster 1 - T. fasciculata") +
  ylab("Median-centered log(CPM)") +
  xlab("Time")

p3b <- ggplot(tfas_cluster1, aes(x=time, y=count, group = gene_id)) +
  geom_line(data = medians_Tfas_cluster1, aes(group = 1), size = 1, color = "black") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  ylim(c(-1.25,-0.25))

p3 <- p3a + annotation_custom(ggplotGrob(p3b), xmin = 5, xmax = 6.5,
                              ymin = -11, ymax = -8)

p4a <- ggplot(tlei_cluster1, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = tlei_cluster1[tlei_cluster1$gene_id == "Tfasc_v1.09028-RA",], aes(group = gene_id),
            size = .8, color = "darkgreen") +
  geom_line(data = tlei_cluster1[tlei_cluster1$gene_id == "Tfasc_v1.19779-RA",], aes(group = gene_id),
            size = .8, color = "darkgreen") +
  geom_line(data = tlei_cluster1[tlei_cluster1$gene_id == "Tfasc_v1.08739-RA",], aes(group = gene_id),
            size = .8, color = "darkgreen") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(ymin_cluster1, ymax_cluster1) +
  ggtitle("2B) cluster 1 - T. leiboldiana") +
  ylab("Median-centered log(CPM)") +
  xlab("Time")

p4b <- ggplot(tlei_cluster1, aes(x=time, y=count, group = gene_id)) +
  geom_line(data = medians_Tlei_cluster1, aes(group = 1), size = 1, color = "black") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

p4 <- p4a + annotation_custom(ggplotGrob(p4b), xmin = 5, xmax = 6.5,
                              ymin = -11, ymax = -8)

######################### CLUSTER 3 #############################
p5a <- ggplot(tfas_cluster3, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = tlei_cluster3[tlei_cluster3$gene_id == "Tfasc_v1.13975-RA",], aes(group = gene_id),
            size = .8, color = "firebrick3") +
  geom_line(data = tlei_cluster3[tlei_cluster3$gene_id == "Tfasc_v1.14652-RA",], aes(group = gene_id),
            size = .8, color = "firebrick3") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(c(ymin_cluster3, ymax_cluster3)) +
  ggtitle("3A) cluster 3 - T. fasciculata") +
  ylab("Median-centered log(CPM)") +
  xlab("Time")

p5b <- ggplot(tfas_cluster3, aes(x=time, y=count, group = gene_id)) +
  geom_line(data = medians_Tfas_cluster3, aes(group = 1), size = 1, color = "black") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

p5 <- p5a + annotation_custom(ggplotGrob(p5b), xmin = 5, xmax = 6.5,
                              ymin = -6.75, ymax = -4.25)

p6a <- ggplot(tlei_cluster3, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = tlei_cluster3[tlei_cluster3$gene_id == "Tfasc_v1.13975-RA",], aes(group = gene_id),
            size = .8, color = "firebrick3") +
  geom_line(data = tlei_cluster3[tlei_cluster3$gene_id == "Tfasc_v1.14652-RA",], aes(group = gene_id),
            size = .8, color = "firebrick3") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ylim(ymin_cluster3, ymax_cluster3) +
  ggtitle("3B) cluster 3 - T. leiboldiana") +
  ylab("Median-centered log(CPM)") +
  xlab("Time")

p6b <- ggplot(tlei_cluster3, aes(x=time, y=count, group = gene_id)) +
  geom_line(data = medians_Tlei_cluster3, aes(group = 1), size = 1, color = "black") +
  geom_vline(xintercept = 4.75, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

p6 <- p6a + annotation_custom(ggplotGrob(p6b), xmin = 5, xmax = 6.5,
                              ymin = -6.75, ymax = -4.25)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2)
dev.off()
