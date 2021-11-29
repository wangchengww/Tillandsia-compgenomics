setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
library(ggplot2)

maptest <- read.delim('mapping_stats.txt', header = T, sep = " ")
summary(maptest)

maptest$ref_genome <- rep(c("A.comosus", "T.fasciculata", "T.leiboldiana"), times = 24)
maptest$sample_species <- c(rep(c("Tfas"), times = 36), rep(c("Tlei"), times = 36))

# Barely any reference bias in T. fasciculata
mean(maptest[maptest$ref_genome == "T.fasciculata" & maptest$sample_species == "Tfas", 3]) # 87.7 %
mean(maptest[maptest$ref_genome == "T.leiboldiana" & maptest$sample_species == "Tfas", 3]) # 87.3 %

# Clear reference bias in T. leiboldiana
mean(maptest[maptest$ref_genome == "T.fasciculata" & maptest$sample_species == "Tlei", 3]) # 84.8 %
mean(maptest[maptest$ref_genome == "T.leiboldiana" & maptest$sample_species == "Tlei", 3]) # 92.9 %

# Multi-mappers: Higher in T. fasciculata than in T. leiboldiana, however, opposite species 
# has higher multimappers to the reference genome.
mean(maptest[maptest$ref_genome == "T.fasciculata",5]) # 7.7 %
mean(maptest[maptest$ref_genome == "T.fasciculata" & maptest$sample_species == "Tlei", 5]) # 8.12 %
mean(maptest[maptest$ref_genome == "T.fasciculata" & maptest$sample_species == "Tfas", 5]) # 7.24 %
mean(maptest[maptest$ref_genome == "T.leiboldiana",5]) # 4.77 %
mean(maptest[maptest$ref_genome == "T.leiboldiana" & maptest$sample_species == "Tlei", 5]) # 4.1 %
mean(maptest[maptest$ref_genome == "T.leiboldiana" & maptest$sample_species == "Tfas", 5]) # 5.48 %


ggplot(data = maptest[maptest$ref_genome == "T.fasciculata",], aes(x=Name, y=uniq, fill=sample_species)) +
  geom_bar(stat="identity") + ylim(c(0,100))
ggplot(data = maptest[maptest$ref_genome == "T.leiboldiana",], aes(x=Name, y=uniq, fill=sample_species)) +
  geom_bar(stat="identity") + ylim(c(0,100))
ggplot(data = maptest[maptest$ref_genome == "A.comosus",], aes(x=Name, y=uniq, fill=sample_species)) +
  geom_bar(stat="identity") + ylim(c(0,100))
