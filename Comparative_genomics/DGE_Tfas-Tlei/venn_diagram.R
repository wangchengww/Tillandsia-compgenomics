setwd("bio-info_phd/DGE_Tfas-Tlei/")

# install.packages("gplots")
library("gplots")

Tfas <- read.table("DE_genes_logFC_2_Tfas_PM-AM.txt", header=F, row.names=NULL)
Tlei <- read.table("DE_genes_logFC_2_Tlei_PM-AM.txt", header=F, row.names=NULL)
AM <- read.table("DE_genes_logFC_2_Tfas-Tlei_AM.txt", header=F, row.names=NULL)
PM <- read.table("DE_genes_logFC_2_Tfas-Tlei_PM.txt", header=F, row.names=NULL)

venn.lists <- list(Tfas, Tlei, AM, PM)
names(venn.lists) <- c("DE_Tfas", "DE_Tlei", "DE_AM", "DE_PM")
v.tableAll <- venn(venn.lists, showSetLogicLabel = F)

venn.lists <- list(Tfas, Tlei)
names(venn.lists) <- c("DE_Tfas", "DE_Tlei")
v.table-species <- venn(venn.lists, showSetLogicLabel = F)

venn.lists <- list(AM, PM)
names(venn.lists) <- c("DE_AM", "DE_PM")
v.table-time <- venn(venn.lists, showSetLogicLabel = F)

# For LogFC 4

Tfas <- read.table("DE_genes_logFC_4_Tfas_PM-AM.txt", header=F, row.names=NULL)
Tlei <- read.table("DE_genes_logFC_4_Tlei_PM-AM.txt", header=F, row.names=NULL)
AM <- read.table("DE_genes_logFC_4_Tfas-Tlei_AM.txt", header=F, row.names=NULL)
PM <- read.table("DE_genes_logFC_4_Tfas-Tlei_PM.txt", header=F, row.names=NULL)

venn.lists <- list(Tfas, Tlei, AM, PM)
names(venn.lists) <- c("DE_Tfas", "DE_Tlei", "DE_AM", "DE_PM")
v.tableAll <- venn(venn.lists, showSetLogicLabel = F)

venn.lists <- list(Tfas, Tlei)
names(venn.lists) <- c("DE_Tfas", "DE_Tlei")
v.table-species <- venn(venn.lists, showSetLogicLabel = F)

venn.lists <- list(AM, PM)
names(venn.lists) <- c("DE_AM", "DE_PM")
v.table-time <- venn(venn.lists, showSetLogicLabel = F)
