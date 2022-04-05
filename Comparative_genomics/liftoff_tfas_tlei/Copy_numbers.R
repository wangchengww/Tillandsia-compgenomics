setwd("/Users/clara/bio-info_phd/liftoff_tfas_tlei/")

d <- read.table('CNV_Tlei_realgenes.txt')
colnames(d) <- c("CN", "geneID")
library("ggplot2")
library(wesanderson)

ggplot(data = d, aes(x=CN)) + geom_histogram()
duplicated <- d[d$CN > 1,]
ggplot(data = duplicated, aes(x=CN)) + 
  geom_histogram(binwidth = 1,color=wes_palette("GrandBudapest1")[3], fill=wes_palette("GrandBudapest1")[4]) +
  labs(x = "Copy number")

ggplot(duplicated, aes(x = CN)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)))

t <- duplicated[duplicated$CN == 2,]
sum(duplicated$CN)
