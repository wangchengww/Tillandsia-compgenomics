library(ggplot2)
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/Supplementary_figures/Assembly_stats/')

tfas_fai <- read.delim("tillandsia_fasciculata_assembly.sorted.fasta.fai", sep = "\t",
                       header = F, colClasses = c(NA, NA, "NULL", "NULL", "NULL"),
                       col.names = c("Scaffold", "length", "NULL", "NULL", "NULL"))

tlei_fai <- read.delim("tillandsia_leiboldiana_assembly.pilon.upper.sorted.fasta.fai", sep = "\t",
                       header = F, colClasses = c(NA, NA, "NULL", "NULL", "NULL"),
                       col.names = c("Scaffold", "length", "NULL", "NULL", "NULL"))

p1a <- ggplot(tfas_fai, 
       aes(x=factor(Scaffold, levels = rev(stats::reorder(Scaffold, length))), 
                     y=length)) + 
  geom_bar(stat="identity", fill="olivedrab") +
  scale_x_discrete(limits=rev) +
  xlab("Assembly scaffolds by decreasing size") +
  ylab("Scaffold length (bp)") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p1b <- ggplot(tfas_fai, 
              aes(x=factor(Scaffold, levels = rev(stats::reorder(Scaffold, length))), 
                  y=length)) + 
  geom_bar(stat="identity", fill="olivedrab") +
  scale_x_discrete(limits=rev(tail(rev(stats::reorder(tfas_fai$Scaffold, tfas_fai$length)),50))) +
  xlab("") +
  ylab("") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p1 <- p1a + annotation_custom(ggplotGrob(p1b), xmin = 300, xmax = 2000,
                        ymin = 2000000, ymax = 30000000)

p2a <- ggplot(tlei_fai, 
       aes(x=factor(Scaffold, levels = rev(stats::reorder(Scaffold, length))), 
           y=length)) + 
  geom_bar(stat="identity", fill="darkgreen") +
  scale_x_discrete(limits=rev) +
  xlab("Assembly scaffolds by decreasing size") +
  ylab("Scaffold length (bp)") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2b <- ggplot(tlei_fai, 
                     aes(x=factor(Scaffold, levels = rev(stats::reorder(Scaffold, length))), 
                         y=length)) + 
  geom_bar(stat="identity", fill="darkgreen") +
  scale_x_discrete(limits=rev(tail(rev(stats::reorder(tfas_fai$Scaffold, tfas_fai$length)),50))) +
  xlab("Assembly scaffolds by decreasing size") +
  ylab("Scaffold length (bp)") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())