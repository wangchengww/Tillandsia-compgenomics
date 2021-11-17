setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/4. Gene_TE_content/")
library(ggplot2)
data <- read.delim("tillandsia_fasciculata_assembly.sorted.fasta.mod.EDTA.TEanno.per_scaffold.abundances.summary", 
                   header = T, sep = '\t')
ggplot(data, aes(x=length, y=Repetitive_content)) + geom_point() + xlim(c(0,1000000))
ggplot(data, aes(x=Repetitive_content)) + geom_density()

mean(data[-(1:25),3]) # On average, small scaffolds have 95 % repetitive content
mean(data[(1:25),3]) # Large scaffolds have on average 65 % repetitive content

data <- data[order(-data$length),]
data$total_bp <- data$length*(0.01*data$Repetitive_content)
l = sum(data[1:25,15])
a = sum(data$total_bp)
perc = l/a
1- perc
# 36 % of all masked basepairs are in smaller contigs

