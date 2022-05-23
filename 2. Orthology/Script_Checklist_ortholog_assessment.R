setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/2. Orthology/")
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/2. Orthology/")


data <- read.delim("checklist_curated_orthologs_Tfas-Tlei.txt", header = T, sep = "\t")
# Count number of genes in T. fasciculata with the same number of exons and the same 
# length as the largest gene in the orthogroup
dim(data[data$diff_largest_exon_count == 0 & 
           startsWith(data$Gene_id, "Tfas"),]) # 19738 (75 %)
dim(data[data$diff_longest_CDS == 1 & 
           startsWith(data$Gene_id, "Tfas"),]) # 18416 (69.9 % %)
# Count number of genes in T. leiboldiana with the same number of exons as the largest
# gene in the orthogroup
dim(data[data$diff_largest_exon_count == 0 & 
           startsWith(data$Gene_id, "Tlei"),]) # 19786 (84 %)
dim(data[data$diff_longest_CDS == 1 & 
           startsWith(data$Gene_id, "Tlei"),]) # 18353 (77.8 % %)

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
tfas_data <- data[startsWith(data$Gene_id, "Tfas"),] 
tlei_data <- data[startsWith(data$Gene_id, "Tlei"),]
# How many multicopy genes are complete in exon count and length in Tfas
tfas_multicopy <-  tfas_data %>% group_by(orthogroup) %>% filter(n()>1) # 12014 genes
length(unique(tfas_multicopy$orthogroup)) # 3331 orthogroups, meaning that there are 8683 secondary copies
dim(tfas_multicopy[tfas_multicopy$diff_largest_exon_count == 0,]) # 5427 (45 %) of which 2096 (24 %) are a secondary copy
dim(tfas_multicopy[tfas_multicopy$diff_longest_CDS == 1,]) # 4105 (34.2%) of which 774 are a secondary copy

# How many multicopy genes are complete in exon count and length in Tlei
tlei_multicopy <-  tlei_data %>% group_by(orthogroup) %>% filter(n()>1) # 8011 genes
length(unique(tlei_multicopy$orthogroup)) # 2311 orthogroups, meaning that 5700 genes are secondary copies
dim(tlei_multicopy[tlei_multicopy$diff_largest_exon_count == 0,]) # 4213 (52.6 %) of which 1902 (34 %) are a secondary copy
dim(tlei_multicopy[tlei_multicopy$diff_longest_CDS == 1,]) # 2780 (34.7 %) of which 469 are a secondary copy

# How many multicopy Tfas genes are expressed in both species, and how does this compare to completeness
dim(tfas_multicopy[tfas_multicopy$expressed_Tfas == "True",]) # 9994 of which 6663 are a secondary copy
dim(tfas_multicopy[tfas_multicopy$expressed_Tfas == "True" & tfas_multicopy$diff_largest_exon_count == 0,]) # 4923 (90 % of complete genes are expressed in Tfas) 
dim(tfas_multicopy[tfas_multicopy$expressed_Tlei == "True" & tfas_multicopy$diff_largest_exon_count == 0,]) # 4844 (89 % of complete genes are expressed in Tlei)

# How many multicopy Tlei genes are expressed in both species
dim(tlei_multicopy[tlei_multicopy$expressed_Tfas == "True",]) # 6194 of which 3883 are a secondary copy
dim(tlei_multicopy[tlei_multicopy$expressed_Tlei == "True",]) # 6443 of which 4112 are a secondary copy

# Make histogram of exon count and length differences per family
p1 <-ggplot(data = tfas_multicopy, aes(x = diff_longest_CDS)) +
  geom_histogram()
p2 <- ggplot(data = tlei_multicopy, aes(x = diff_longest_CDS)) +
  geom_histogram()
grid.arrange(p1, p2, nrow = 1)

p1 <-ggplot(data = tfas_multicopy, aes(x = diff_largest_exon_count)) +
  geom_histogram()
p2 <- ggplot(data = tlei_multicopy, aes(x = diff_largest_exon_count)) +
  geom_histogram()
grid.arrange(p1, p2, nrow = 1)

# Make same histogram only for expressed genes
p1 <-ggplot(data = tfas_multicopy[tfas_multicopy$expressed_Tfas == "True",], aes(x = diff_longest_CDS)) +
  geom_histogram()
p2 <- ggplot(data = tlei_multicopy[tlei_multicopy$expressed_Tlei == "True",], aes(x = diff_longest_CDS)) +
  geom_histogram()
grid.arrange(p1, p2, nrow = 1)

# How many multicopy genes have both a start & stopcodon, either or none in Tfas
dim(tfas_multicopy[tfas_multicopy$startcodon == "True" & tfas_multicopy$stopcodon == "True",]) # 10,331 (86 % of multi copy genes)
dim(tfas_multicopy[(tfas_multicopy$startcodon == "True" & tfas_multicopy$stopcodon == "False") | 
                     (tfas_multicopy$startcodon == "False" & tfas_multicopy$stopcodon == "True"),]) # 1373 (11.4 % of multicopy genes)
dim(tfas_multicopy[tfas_multicopy$startcodon == "False" & tfas_multicopy$stopcodon == "False",]) # 310 genes (2.58 % of multicopy genes)

# How many multicopy genes have both a start & stopcodon, either or none in Tlei
dim(tlei_multicopy[tlei_multicopy$startcodon == "True" & tlei_multicopy$stopcodon == "True",]) # 4011 (50 % of multicopy genes)
dim(tlei_multicopy[(tlei_multicopy$startcodon == "True" & tlei_multicopy$stopcodon == "False") | 
                     (tlei_multicopy$startcodon == "False" & tlei_multicopy$stopcodon == "True"),]) # 3382 (42 %)
dim(tlei_multicopy[tlei_multicopy$startcodon == "False" & tlei_multicopy$stopcodon == "False",]) # 618 (7.71 %)
