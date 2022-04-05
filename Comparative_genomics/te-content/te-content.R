setwd('/Users/clara/bio-info_phd/Comparative_genomics/te-content/')

te <- read.table("summary_TE_GCcontent_Leiboldiana.txt",header = F, sep = "\t")
colnames(te) <- c("gene_id", "i", "GC", "GC-TE", "ATGC", "ATGC-TE", "N", "GC_perc",
                  "N_perc", "TE_perc", "GC-TE_perc")

s <- te[,c(1,10)]
write.table(s, "list_te_content_per_gene.txt")

te_50 <- te[te$TE_perc >= 0.5,]
te_75 <- te[te$TE_perc >= 0.75,]
te_25 <- te[te$TE_perc >= 0.25,]
te_10 <- te[te$TE_perc <= 0.1,]
te_0 <- te[te$TE_perc <= 0,]
library(ggplot2)
ggplot(te, aes(TE_perc))  +
  geom_density() +
  theme_bw()+ 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_vline(xintercept=0.75,col="red")+
  geom_vline(xintercept=0.5,col="orange")+
  geom_vline(xintercept=0.25,col="darkgreen")+
  annotate("text", x=0.84, y=5, label= "1947 genes",col="red")+
  annotate("text", x=0.60, y=5, label= "3896 genes",col="orange")+
  annotate("text", x=0.35, y=5, label= "7067 genes",col="darkgreen")+
  xlab("Masked percentage of CDS")+
  ylab("Gene density")+
  labs(title = "Te content of T. leiboldiana gene models (maker R1)")
