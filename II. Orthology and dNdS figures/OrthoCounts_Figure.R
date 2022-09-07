# TL - 050922
# CGC - 070922
library(ggplot2)
library(dplyr)
library(hexbin)
library(ggrepel)
library(grid)
library(gridExtra)

setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/II. Orthology and dNdS figures/")
counts <- read.table("orthogroups_Tfas_Tlei_Acom.counts.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.blastandsearch.txt", sep = '\t')
colnames(counts) <- c("og_id", "Acom", "Tfas", "Tlei")

# Filter out unique orthogroups
counts_Tfas_Tlei <- counts[counts$Tfas != 0 & counts$Tlei != 0,]
# Remove 1:1 orthogroups
counts_Tfas_Tlei_multi <- counts_Tfas_Tlei[!(counts_Tfas_Tlei$Tfas == 1 & counts_Tfas_Tlei$Tlei == 1),]

### Make figure ###
hist_top <- ggplot() + 
  geom_histogram(color="black", fill="olivedrab",
                 aes(counts_Tfas_Tlei_multi$Tfas),
                 bins=ceiling(max(counts_Tfas_Tlei_multi$Tfas)))+
  xlab(label = "T. fasciculata") +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="olivedrab",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="plain"))

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

scatter <- ggplot(counts_Tfas_Tlei_multi, aes(x=Tfas, y=Tlei)) + geom_point(size = .5) +
  geom_abline(intercept = 0, slope = 1) +
  ylab(label = "T. leiboldiana") +
  xlab(label = "T. fasciculata")+
  scale_fill_manual(values=c("darkgreen","olivedrab")) +
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="olivedrab",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="darkgreen",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))+
  geom_label_repel(data=subset(counts_Tfas_Tlei_multi, 
                               abs(counts_Tfas_Tlei_multi$Tfas - counts_Tfas_Tlei_multi$Tlei)>=15),
                   aes(label = og_id,
                       fill=ifelse(Tlei>=Tfas,"darkgreen","grey40")),
                   colour = 'white',size = 3)

hist_right <- ggplot()+ 
  geom_histogram(color="black", fill="darkgreen",
                 aes(counts_Tfas_Tlei_multi$Tlei),
                bins=ceiling(max(counts_Tfas_Tlei_multi$Tlei)))+
  coord_flip() +
  xlab(label = "T. leiboldiana") +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="plain"),
        axis.title.y = element_text(colour="darkgreen",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


# print results !
pdf("Figure_OrthogroupCounts_GeneFamilyEvo.pdf", height = 8, width = 8)
grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, 
             widths=c(4, 1),heights=c(1, 4))
dev.off()
