setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
library(ggplot2)
library(grid)
library(gridExtra)

maptest <- read.delim('mapping_stats.testset.txt', header = T, sep = " ")
mask <- read.delim('mapping_stats.testset.masked.txt', header = T, sep = ' ')
summary(maptest)
summary(mask)

maptest$ref_genome <- rep(c("A.comosus", "T.fasciculata", "T.leiboldiana"), times = 24)
maptest$sample_species <- c(rep(c("Tfas"), times = 36), rep(c("Tlei"), times = 36))

mask$ref_genome <- rep(c("T.fasciculata", "T.leiboldiana"), times = 24)
mask$sample_species <- c(rep(c("Tfas"), times = 24), rep(c("Tlei"), times = 24))

# Reference bias in T. fasciculata: 2.9 % difference on average
mean(maptest[maptest$ref_genome == "T.fasciculata" & maptest$sample_species == "Tfas", 3]) # 87.7 %
mean(maptest[maptest$ref_genome == "T.fasciculata" & maptest$sample_species == "Tlei", 3]) # 84.8 %

mean(mask[mask$ref_genome == "T.fasciculata" & mask$sample_species == "Tfas", 3]) # 85.8 %
mean(mask[mask$ref_genome == "T.fasciculata" & mask$sample_species == "Tlei", 3]) # 83.1 %

# Reference bias in T. leiboldiana: 5.6 % difference on average
mean(maptest[maptest$ref_genome == "T.leiboldiana" & maptest$sample_species == "Tlei", 3]) # 92.9 %
mean(maptest[maptest$ref_genome == "T.leiboldiana" & maptest$sample_species == "Tfas", 3]) # 87.3 %

mean(mask[mask$ref_genome == "T.leiboldiana" & mask$sample_species == "Tlei", 3]) # 90.5 %
mean(mask[mask$ref_genome == "T.leiboldiana" & mask$sample_species == "Tfas", 3]) # 85.8 %

# In A.comosus
mean(maptest[maptest$ref_genome == "A.comosus" & maptest$sample_species == "Tlei", 3]) # 39.83 %
mean(maptest[maptest$ref_genome == "A.comosus" & maptest$sample_species == "Tfas", 3]) # 39.9 %

# Multi-mappers: Higher in T. fasciculata than in T. leiboldiana, however, opposite species 
# has higher multimappers to the reference genome.
mean(maptest[maptest$ref_genome == "T.fasciculata",5]) # 7.7 %
mean(maptest[maptest$ref_genome == "T.fasciculata" & maptest$sample_species == "Tlei", 5]) # 8.12 %
mean(maptest[maptest$ref_genome == "T.fasciculata" & maptest$sample_species == "Tfas", 5]) # 7.24 %
mean(maptest[maptest$ref_genome == "T.leiboldiana",5]) # 4.77 %
mean(maptest[maptest$ref_genome == "T.leiboldiana" & maptest$sample_species == "Tlei", 5]) # 4.1 %
mean(maptest[maptest$ref_genome == "T.leiboldiana" & maptest$sample_species == "Tfas", 5]) # 5.48 %


p1 <- ggplot(data = maptest[maptest$ref_genome == "T.fasciculata",], aes(x=Name, y=uniq, fill=sample_species)) +
  geom_bar(stat="identity") + 
  xlab("") +
  scale_y_continuous(name="% reads mapping uniquely", 
                     breaks = c(seq(from = 10, to = 110, by = 10)), 
                     limits = c(0,100)) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "none")

p2 <- ggplot(data = maptest[maptest$ref_genome == "T.leiboldiana",], aes(x=Name, y=uniq, fill=sample_species)) +
  geom_bar(stat="identity") +
  xlab("Samples") + 
  scale_y_continuous(name="", 
                     breaks = c(seq(from = 10, to = 100, by = 10)), 
                     limits = c(0,100)) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "none")

p3 <- ggplot(data = maptest[maptest$ref_genome == "A.comosus",], aes(x=Name, y=uniq, fill=sample_species)) +
  geom_bar(stat="identity") + xlab("") +
  scale_y_continuous(name="", 
                     breaks = c(seq(from = 10, to = 100, by = 10)), 
                     limits = c(0,100)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position = "none")

grid.arrange(p1, p2, p3, nrow = 1)

# For masked alignments
p4 <- ggplot(data = mask[mask$ref_genome == "T.fasciculata",], aes(x=Name, y=uniq, fill=sample_species)) +
  geom_bar(stat="identity") + 
  xlab("") +
  scale_y_continuous(name="% reads mapping uniquely", 
                     breaks = c(seq(from = 10, to = 110, by = 10)), 
                     limits = c(0,100)) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "none")

p5 <- ggplot(data = mask[mask$ref_genome == "T.leiboldiana",], aes(x=Name, y=uniq, fill=sample_species)) +
  geom_bar(stat="identity") +
  xlab("Samples") + 
  scale_y_continuous(name="", 
                     breaks = c(seq(from = 10, to = 100, by = 10)), 
                     limits = c(0,100)) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "none")

grid.arrange(p4, p5, nrow = 1)
grid.arrange(p1, p2, p4, p5, nrow = 1)

#------------------- GENERAL MAPPING STATS ------------------------

map_data <- read.delim("mapping_stats.full.txt", sep = " ", header = T)
mean(map_data$total) # 68,856,922 reads per sample on average
min(map_data$total) # 49,469,111 minimum reads
quantile(map_data$total, c(0.05,0.95)) # 90 % of samples have between 53 - 82 million reads

mean(map_data$uniq) # 85.4 % of reads per sample are uniquely mapped on average
map_data$sample_species <- c(rep(c("Tfas"), times = 36), rep(c("Tlei"), times = 36))
map_data$sample_number <- c(rep(c("A"), times = 6), rep(c("B"), times = 6),
                            rep(c("C"), times = 6), rep(c("D"), times = 6),
                            rep(c("E"), times = 6), rep(c("F"), times = 6),
                            rep(c("A"), times = 6), rep(c("C"), times = 6),
                            rep(c("D"), times = 6), rep(c("E"), times = 6),
                            rep(c("F"), times = 6), rep(c("G"), times = 6))
                    

mean(map_data[map_data$sample_species == "Tfas", 3]) # 85.7 %
mean(map_data[map_data$sample_species == "Tlei", 3]) # 85.19 %
min(map_data$uniq)

p6 <- ggplot(data = map_data, aes(x=Name, y=uniq, fill=sample_number)) +
  geom_bar(stat="identity") + 
  xlab("Samples") +
  scale_y_continuous(name="% reads mapping uniquely", 
                     breaks = c(seq(from = 0, to = 110, by = 5))) +
  coord_cartesian(ylim=c(65,100)) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "none") +
  geom_hline(yintercept=85.4, linetype = "dashed", size = .2)
  

grid.arrange(p4, p5, nrow = 1)



