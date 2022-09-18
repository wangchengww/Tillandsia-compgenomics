pacman::p_load("ggplot2","grid", "gridExtra")
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/Supplementary_figures/Assembly_stats/')

tfas_fai <- read.delim("tillandsia_fasciculata_assembly.sorted.fasta.fai", sep = "\t",
                       header = F, colClasses = c(NA, NA, "NULL", "NULL", "NULL"),
                       col.names = c("Scaffold", "length", "NULL", "NULL", "NULL"))

tlei_fai <- read.delim("tillandsia_leiboldiana_assembly.pilon.upper.sorted.fasta.fai", sep = "\t",
                       header = F, colClasses = c(NA, NA, "NULL", "NULL", "NULL"),
                       col.names = c("Scaffold", "length", "NULL", "NULL", "NULL"))

# Separate plots for Tfas and Tlei
p1 <- ggplot(tfas_fai, 
       aes(x=factor(Scaffold, levels = rev(stats::reorder(Scaffold, length))), 
                     y=length)) + 
  geom_bar(stat="identity", fill="olivedrab") +
  scale_x_discrete(limits=rev(tail(rev(stats::reorder(tfas_fai$Scaffold, tfas_fai$length)),100))) +
  xlab("First 100 scaffolds by decreasing size") +
  ylim(c(0,60000000)) +
  ylab("Scaffold length (bp)") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2 <- ggplot(tlei_fai, 
       aes(x=factor(Scaffold, levels = rev(stats::reorder(Scaffold, length))), 
           y=length)) + 
  geom_bar(stat="identity", fill="darkgreen") +
  scale_x_discrete(limits=rev(tail(rev(stats::reorder(tlei_fai$Scaffold, tlei_fai$length)),100))) +
  xlab("First 100 scaffolds by decreasing size") +
  ylab("Scaffold length (bp)") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

grid.arrange(p1, p2, nrow = 2, ncol = 1)

# Joined plot
tfas_fai$species <- 'tfas'
tfas_fai$order <- order(tfas_fai$length, decreasing = T)
tlei_fai$species <- 'tlei'
tlei_fai$order <- order(tlei_fai$length, decreasing = T)
all_fai <- rbind(tlei_fai, tfas_fai)
pdf("SuppFigure_AssemblyStats_ScaffoldSizes.pdf", width = 9, height = 6)
p3 <- ggplot(all_fai, 
             aes(x=order, y=length, fill= species)) + 
  geom_bar(stat="identity", position = "identity", alpha = .7) +
  xlim(c(0,100)) +
  theme_bw() +
  xlab("First 100 scaffolds by decreasing size") +
  ylab("Scaffold length (bp)") +
  scale_fill_manual(values=c("grey40", "darkgreen"), name = "Species", 
                      labels = c("T. fasciculata", "T. leiboldiana")) +
  theme(legend.text = element_text(face = "italic"))
p3
dev.off()

# Genes per scaffold

count_table_Tfas <- read.delim("Tfas_perScaffold_Orthogenes_count.txt", sep = "\t",
                               header = T)
count_table_Tlei <- read.delim("Tlei_perScaffold_Orthogenes_count.txt", sep = "\t",
                               header = T)
count_table_Tfas$species <- 'tfas'
count_table_Tlei$species <- 'tlei'
all_genes <- rbind(count_table_Tfas, count_table_Tlei)

pdf("SuppFigure_AssemblyStats_GeneCounts_perScaff.pdf", width = 9, height = 6)
p4 <- ggplot(all_genes, 
             aes(x=order, y=orthocount_all, fill= species)) + 
  geom_bar(stat="identity", position = "identity", alpha = .7) +
  xlim(c(0,100)) +
  theme_bw() +
  xlab("First 100 scaffolds by decreasing size") +
  ylab("Number of orthologous genes") +
  scale_fill_manual(values=c("grey40", "darkgreen"), name = "Species", 
                    labels = c("T. fasciculata", "T. leiboldiana")) +
  theme(legend.text = element_text(face = "italic"))
p4
dev.off()

pdf("SuppFigure_AssemblyStats_COMBINED.pdf", width = 12, height = 6)
grid.arrange(p3, p4, nrow = 2, ncol = 1)
dev.off()
