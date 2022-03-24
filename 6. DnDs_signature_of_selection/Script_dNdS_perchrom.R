setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/6. DnDs_signature_of_selection/')

pacman::p_load("ggplot2", "gtools","grid", "gridExtra")

data <- read.table("pairwise_dnds.codeml.allgenes.MACSE.with-chrom.output", sep = "\t", header = T)
summary(data)
data$chrom_Tfas <- factor(data$chrom_Tfas, levels=mixedsort(unique(data$chrom_Tfas)))
data$chrom_Tlei <- factor(data$chrom_Tlei, levels=mixedsort(unique(data$chrom_Tlei)))

data_sign <- data[data$pvalue < 0.05,]
mean(data_sign$dN.dS) 

mean(data$dN.dS) # 1.986
median(data$dN.dS) # 0.347
mean(data$dN)
mean(data$dS)

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
  geom_boxplot() + ylim(c(0,1))

