#clean masigpro history
#Karolina Heyduk, University of Georgia, 2016
#modified 2022 (annotations)
#contact: heyduk@hawaii.edu

#load data and libraries
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/")
counts <- read.table("../counts.Tfas.6_timepoints.txt", header = T, row.names = 1) #important that these are the raw counts
library(edgeR)
library(maSigPro)
library(mclust)
library(stringr)

#set up edgeR object
groups <- unlist(data.table::transpose(str_split(colnames(counts), "_"))[c(4)])
dyg<-DGEList(counts, group=groups)
dyg<-calcNormFactors(dyg, method="TMM")
normd<-cpm(dyg, normalize.lib.sizes=TRUE)
plotMDS(normd, labels=groups) #lets you look for outliers and remove if necessary

##remove genes that have all zero read counts
normd<-normd[ rowSums(normd)!=0, ]

##load design object for masigpro. row order must be the same as the order of libraries in count matrix
design <- data.frame(Time=c(rep(c(1,2,3,4,5,6))), 
                     Replicates=c(rep(1,6), rep(2,6), rep(3,6), rep(4,6),rep(5,6), rep(6,6)),
                     Tfas = c(rep(1, 36)))
rownames(design) <- colnames(counts)
design
d<-make.design.matrix(design, degree=5)
##using the negative binomial options in masigpro, calculate polynomial regressions for each gene
NBp<-p.vector(normd, d, counts=TRUE) #please choose dis. family with care. default for counts is neg. binomial
NBp$i # returns the number of significant genes

##TO REMOVE INFLUENTIAL GENES:
NBt<-T.fit(NBp)
influential<-NBt$influ.info
inf.genenames<-colnames(influential)
normd<-normd[!rownames(normd) %in% inf.genenames, ]

# Get significant genes 
sigs <- get.siggenes(NBt, rsq = 0.6, vars = "groups")
sigs$summary
write.table(sigs$summary, file = "Genes_Significant_Tfas.txt", quote = F, sep = "\t", row.names = F)
see.genes(sigs$sig.genes$Tfas, show.fit = T, dis=d$dis,
          cluster.method="Mclust", cluster.data = 1, k.mclust = T)

  ##pick k
NBp<-p.vector(normd, d, counts=TRUE)
wss<-(nrow(NBp$SELEC)-1)*sum(apply(NBp$SELEC,2,var))
for (i in 2:15) wss[i]<- sum(kmeans(NBp$SELEC, centers=i, iter.max=20)$withinss)
plot(1:15, wss, type="b")


