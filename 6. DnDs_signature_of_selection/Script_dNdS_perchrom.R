setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/6. DnDs_signature_of_selection/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/6. DnDs_signature_of_selection/')

pacman::p_load("ggplot2", "gtools","grid", "gridExtra", "dplyr")

data <- read.table("pairwise_dnds.codeml.allgenes.MACSE.with-chrom.filtered-5variants.output", sep = "\t", header = T)
summary(data)
data$chrom_Tfas <- factor(data$chrom_Tfas, levels=mixedsort(unique(data$chrom_Tfas)))
data$chrom_Tlei <- factor(data$chrom_Tlei, levels=mixedsort(unique(data$chrom_Tlei)))

data_sign <- data[data$pvalue < 0.05,]
mean(data_sign$dN.dS) 

mean(data$dN.dS) # 1.986
median(data$dN.dS) 
mean(data$dN)
mean(data$dS)

avgs_Tfas <- data %>%
  group_by(chrom_Tfas) %>%
  dplyr::summarize(Median = median(dN.dS, na.rm=TRUE))

avgs_Tlei <- data %>%
  group_by(chrom_Tlei) %>%
  dplyr::summarize(Median = median(dN.dS, na.rm=TRUE))

ggplot(data, aes(x=dN.dS)) + geom_density() +
  xlim(c(0,6)) + theme_linedraw()

ggplot(data, aes(x=dN.dS)) + geom_histogram() +
  xlim(c(0,6)) + theme_linedraw()

p1 <- ggplot(data, aes(x=(chrom_Tfas), y=dN.dS)) + 
  geom_boxplot() + ylim(c(0,2))
p2 <- ggplot(data, aes(x=(chrom_Tlei), y=dN.dS)) + 
  geom_boxplot() + ylim(c(0,2))

grid.arrange(p1, p2, nrow = 2, top = "Per-chromosome dN/dS ratio distribution for T. fasciculata and T. leiboldiana")

ggplot(data_sign, aes(x=(chrom_Tfas), y=dN.dS)) + 
  geom_boxplot() + ylim(c(0,2))

# Plot dNdS values across a chromosome with a fusion (chr14)
# The breakpoint was determined to be around 28,588,630 – 28,625,508 bp.
data_chr14 <- data[data$chrom_Tlei == 'chr14',]
data_sign_chr14 <- data_sign[data_sign$chrom_Tlei == 'chr14',]
pdf("Tlei-chr14_dNdS_across_chromosome-strict_ZOOM.pdf", width = 8, height = 6)
ggplot(data_sign_chr14, aes(x=start.1, y=dN.dS)) + geom_point() +
  ylim(c(0,5)) +
  geom_rect(aes(xmin=28588630, xmax=28625508, ymin=0, ymax=5), color="red") +
  ggtitle("dNdS values across chromosome 14 in T. leiboldiana\nZOOMED-IN\nonly significant genes in MACSE") +
  ylab("dN/dS") +
  xlab("position in chr14")
dev.off()
# Now for Tfas chromosome 10 where a translocation has occurred
# Breakpoint: 22,956,800 – 22,963,481
data_chr10 <- data[data$chrom_Tfas == 'chr10',]
data_sign_chr10 <- data_sign[data_sign$chrom_Tfas == 'chr10',]
ggplot(data_sign_chr10, aes(x=start, y=dN.dS)) + geom_point() +
  geom_rect(aes(xmin=22956800, xmax=22963481, ymin=0, ymax=5), color="red")
