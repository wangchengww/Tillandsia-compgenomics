##
# Set working directory and read in data
setwd('/Users/clara/bio-info_phd/Comparative_genomics/orthology_Tfas_Tlei/synteny_all_scaffolds/reanalysis_august2021_newcounts/')
ortho_tfas <- read.delim("Tfas_orthology_info.per_scaffold.txt", header = F, sep = "\t")
colnames(ortho_tfas) <- c('geneID', 'scaffold', 'scaffold_length', 'startpos','endpos','ogID','Acom_count', 'Tfas_count', 'Tlei_count', 'GOterms', 'Description')
# Calculate the length of each gene
ortho_tfas$length_gene <- ortho_tfas$endpos - ortho_tfas$startpos
# Same for Tlei
ortho_tlei <- read.delim("Tlei_orthology_info.per_scaffold.txt", header = F, sep = "\t")
colnames(ortho_tlei) <- c('geneID', 'scaffold', 'scaffold_length', 'startpos','endpos','ogID','Acom_count', 'Tfas_count', 'Tlei_count', 'GOterms', 'Description')
ortho_tlei$length_gene <- ortho_tlei$endpos - ortho_tlei$startpos

##
# Calculate the number of orthologous genes on each scaffold (all relations)
count_per_scaffold_Tfas <- data.frame(matrix(ncol=3))
c <- c()
for (scaffold in unique(ortho_tfas$scaffold)){
  count <- (nrow(ortho_tfas[ortho_tfas$scaffold == paste0(scaffold),]))
  length <- unique(ortho_tfas[ortho_tfas$scaffold == paste0(scaffold),3])
  c <- c(scaffold,length, count)
  count_per_scaffold_Tfas <- rbind(count_per_scaffold_Tfas,c)
  colnames(count_per_scaffold_Tfas) <- c("scaffold", "length", "ortholog_count")
}
count_per_scaffold_Tfas <- count_per_scaffold_Tfas[-1,]
# Same for Tlei
count_per_scaffold_Tlei <- data.frame(matrix(ncol=3))
c <- c()
for (scaffold in unique(ortho_tlei$scaffold)){
  count <- (nrow(ortho_tlei[ortho_tlei$scaffold == paste0(scaffold),]))
  length <- unique(ortho_tlei[ortho_tlei$scaffold == paste0(scaffold),3])
  c <- c(scaffold,length, count)
  count_per_scaffold_Tlei <- rbind(count_per_scaffold_Tlei,c)
  colnames(count_per_scaffold_Tlei) <- c("scaffold", "length", "ortholog_count")
}
count_per_scaffold_Tlei <- count_per_scaffold_Tlei[-1,]

##
# Select orthologues with 1:1 relation
oneone_Tfas <- ortho_tfas[(ortho_tfas$Tfas_count==1 & (ortho_tfas$Tlei_count==1)),]
# Calculate count of 1:1 orthologs per scaffold
oneone_count_per_scaffold_Tfas <- data.frame(matrix(ncol=3))
c <- c()
for (scaffold in unique(oneone_Tfas$scaffold)){
  count <- (nrow(oneone_Tfas[oneone_Tfas$scaffold == paste0(scaffold),]))
  length <- unique(oneone_Tfas[oneone_Tfas$scaffold == paste0(scaffold),3])
  c <- c(scaffold,length, count)
  oneone_count_per_scaffold_Tfas <- rbind(oneone_count_per_scaffold_Tfas,c)
  colnames(oneone_count_per_scaffold_Tfas) <- c("scaffold", "length", "ortholog_count")
}
oneone_count_per_scaffold_Tfas <- oneone_count_per_scaffold_Tfas[-1,]
#Same for Tlei
oneone_Tlei <- ortho_tlei[(ortho_tlei$Tfas_count==1 & (ortho_tlei$Tlei_count==1)),]
oneone_count_per_scaffold_Tlei <- data.frame(matrix(ncol=3))
c <- c()
for (scaffold in unique(oneone_Tlei$scaffold)){
  count <- (nrow(oneone_Tlei[oneone_Tlei$scaffold == paste0(scaffold),]))
  length <- unique(oneone_Tlei[oneone_Tlei$scaffold == paste0(scaffold),3])
  c <- c(scaffold,length, count)
  oneone_count_per_scaffold_Tlei <- rbind(oneone_count_per_scaffold_Tlei,c)
  colnames(oneone_count_per_scaffold_Tlei) <- c("scaffold", "length", "ortholog_count")
}
oneone_count_per_scaffold_Tlei <- oneone_count_per_scaffold_Tlei[-1,]

##
# Merge 1:1 counts with all orthologous counts and compute ratio between all counts and 1:1 counts
c2 <- merge(count_per_scaffold_Tfas, oneone_count_per_scaffold_Tfas, by = "scaffold", all = T)
c2 <- c2[,c(1:3,5)]
c2[is.na(c2)] <- 0
count_table_Tfas <- c2
colnames(count_table_Tfas) <- c("scaffold", "length", "orthocount_all", "orthocount_one")
count_table_Tfas$length <- as.numeric(count_table_Tfas$length)
count_table_Tfas$orthocount_all <- as.numeric(count_table_Tfas$orthocount_all)
count_table_Tfas$orthocount_one <- as.numeric(count_table_Tfas$orthocount_one)
count_table_Tfas$ratio <- (count_table_Tfas$orthocount_one)/(count_table_Tfas$orthocount_all)
count_table_Tfas <- count_table_Tfas[order(count_table_Tfas$length, decreasing = T),]
row.names(count_table_Tfas) <- NULL #resets row index
count_table_Tfas$order <- c(1:1155)
# Plot the number of orthologous genes per scaffold
# We use tidyr to transform the dataset so that the "orthocount" columns become
# observations in the new column "type". The corresponding count is placed in a column
# called "count"
library(ggplot2)
library(tidyr)
count_table_Tfas[1:30,] %>%
  gather(type, count, orthocount_all:orthocount_one) %>%
  ggplot(., aes(x=order, y=count, fill=forcats::fct_rev(type))) +
  geom_bar(stat="identity", position="identity") +
  xlab("Scaffolds in decreasing order by size") +
  ylab("Number of orthologous genes") +
  labs(title = "Distribution of orthologous genes in T. fasciculata assembly") +
  scale_fill_discrete(name = "Orthology", labels = c("1:1", "all"))

# Same for Tlei
c2lei <- merge(count_per_scaffold_Tlei, oneone_count_per_scaffold_Tlei, by = "scaffold", all = T)
c2lei <- c2lei[,c(1:3,5)]
c2lei[is.na(c2lei)] <- 0
count_table_Tlei <- c2lei
colnames(count_table_Tlei) <- c("scaffold", "length", "orthocount_all", "orthocount_one")
count_table_Tlei$length <- as.numeric(count_table_Tlei$length)
count_table_Tlei$orthocount_all <- as.numeric(count_table_Tlei$orthocount_all)
count_table_Tlei$orthocount_one <- as.numeric(count_table_Tlei$orthocount_one)
count_table_Tlei$ratio <- (count_table_Tlei$orthocount_one)/(count_table_Tlei$orthocount_all)
count_table_Tlei <- count_table_Tlei[order(count_table_Tlei$length, decreasing = T),]
row.names(count_table_Tlei) <- NULL #resets row index
count_table_Tlei$order <- c(1:2210)
# Plot the number of orthologous genes per scaffold
library(ggplot2)
library(tidyr)
count_table_Tlei[1:30,] %>%
  gather(type, count, orthocount_all:orthocount_one) %>%
  ggplot(., aes(x=order, y=count, fill=forcats::fct_rev(type))) +
  geom_bar(stat="identity", position="identity") +
  xlab("Scaffolds in decreasing order by size") +
  ylab("Number of orthologous genes") +
  labs(title = "Distribution of orthologous genes in T. leiboldiana assembly") +
  scale_fill_discrete(name = "Orthology", labels = c("1:1", "all"))

##
# Number of orthologous genes in 25 largest scaffolds
subset_tfas <- count_table_Tfas[1:25,]
sum(subset_tfas$orthocount_all) # 28497 (90.3 %)
sum(count_table_Tfas$orthocount_all) # 31551
sum(subset_tfas$orthocount_one) # 12667 (99 %)
sum(count_table_Tfas$orthocount_one) # 12783
# Same for Tlei
subset_tlei <- count_table_Tlei[1:25,]
sum(subset_tlei$orthocount_all) # 28859 (88.6 %)
sum(count_table_Tlei$orthocount_all) # 32568
sum(subset_tlei$orthocount_one) # 12725 (99.5 %)
sum(count_table_Tlei$orthocount_one) # 12783

##
# Length distribution
subset_count_Tlei <- ortho_tlei[(ortho_tlei$Tfas_count==1 & (ortho_tlei$Tlei_count==1)),]
subset_count_Tlei <- rbind(subset_count_Tlei, ortho_tlei[(ortho_tlei$Tfas_count==0 & (ortho_tlei$Tlei_count==1)),])
subset_count_Tlei <- rbind(subset_count_Tlei, ortho_tlei[(ortho_tlei$Tfas_count==1 & (ortho_tlei$Tlei_count==2)),])
class <- c()
for (i in 1:nrow(subset_count_Tlei)){
  if (subset_count_Tlei[i,7] == 1 & subset_count_Tlei[i,8] == 1){
    x = "1:1"
    class = c(class, x)
  } else if (subset_count_Tlei[i,7] == 0 & subset_count_Tlei[i,8] == 1) {
    x = "0:1"
    class = c(class, x)
  } else if (subset_count_Tlei[i,7] == 1 & subset_count_Tlei[i,8] == 2) {
    x = "1:2"
    class = c(class, x)
  }
}
subset_count_Tlei$class <- class
library(plyr)
mu <- ddply(subset_count_Tlei, "class", summarise, grp.mean=mean(length_gene))
head(mu)
library(wesanderson)
palette_Tlei <- c(wes_palette("Darjeeling1")[1], wes_palette("Darjeeling1")[2],wes_palette("Darjeeling1")[4])
ggplot(subset_count_Tlei, aes(x=length_gene, fill=class)) +
  geom_density(alpha=0.8) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=class),
             linetype="dashed") +
  scale_fill_manual(values=palette_Tlei) +
  scale_color_manual(values=palette_Tlei) +
  labs(title = "Gene length distribution for three classes of orthogroups in T. leiboldiana") +
  xlab(label = "Gene length (bp)")
# Same for Tfas
subset_count_Tfas <- ortho_tfas[(ortho_tfas$Tfas_count==1 & (ortho_tfas$Tlei_count==1)),]
subset_count_Tfas <- rbind(subset_count_Tfas, ortho_tfas[(ortho_tfas$Tfas_count==1 & (ortho_tfas$Tlei_count==0)),])
subset_count_Tfas <- rbind(subset_count_Tfas, ortho_tfas[(ortho_tfas$Tfas_count==2 & (ortho_tfas$Tlei_count==1)),])
class <- c()
for (i in 1:nrow(subset_count_Tfas)){
  if (subset_count_Tfas[i,7] == 1 & subset_count_Tfas[i,8] == 1){
    x = "1:1"
    class = c(class, x)
  } else if (subset_count_Tfas[i,7] == 1 & subset_count_Tfas[i,8] == 0) {
    x = "1:0"
    class = c(class, x)
  } else if (subset_count_Tfas[i,7] == 2 & subset_count_Tfas[i,8] == 1) {
    x = "2:1"
    class = c(class, x)
  }
}
subset_count_Tfas$class <- class
library(plyr)
mu <- ddply(subset_count_Tfas, "class", summarise, grp.mean=mean(length_gene))
head(mu)
library(wesanderson)
palette_Tfas <- c(wes_palette("Darjeeling2")[1], wes_palette("Darjeeling2")[2],wes_palette("Darjeeling2")[4])
ggplot(subset_count_Tfas, aes(x=length_gene, fill=class)) +
  geom_density(alpha=0.8) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=class),
             linetype="dashed") +
  scale_fill_manual(values=palette_Tfas) +
  scale_color_manual(values=palette_Tfas) +
  labs(title = "Gene length distribution for three classes of orthogroups in T. fasciculata") +
  xlab(label = "Gene length (bp)")
