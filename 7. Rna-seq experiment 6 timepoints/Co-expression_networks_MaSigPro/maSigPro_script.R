#clean masigpro history
#Karolina Heyduk, University of Georgia, 2016
#modified 2022 (annotations)
#contact: heyduk@hawaii.edu

#load data and libraries
setwd("C:/Users/heyduk/Google Drive/Science/CAM transcription/2016/maProSig/")
counts<-read.table("../monocot ortho/YG.counts.matrix") #important that these are the raw counts
library(edgeR)
library(maSigPro)
library(Mfuzz)

#set up edgeR object
groups<-c("YG.2.D","YG.6.D","YG.1.D","YG.5.D","YG.4.W","YG.1.D","YG.4.W","YG.2.W","YG.6.W","YG.6.D","YG.4.D","YG.5.D","YG.1.W","YG.5.W","YG.1.W","YG.5.D","YG.3.D","YG.1.W","YG.5.W","YG.3.D","YG.5.D","YG.4.D","YG.5.W","YG.6.W","YG.1.D","YG.6.D","YG.2.W","YG.4.W","YG.2.W","YG.3.W","YG.3.D","YG.3.W","YG.2.D","YG.3.W","YG.4.W","YG.5.W","YG.3.W","YG.1.D","YG.4.D","YG.3.D","YG.2.W","YG.2.D","YG.4.D","YG.2.D","YG.6.D") #this will be specific to your data. Not necessarily required, but nice for checking that replicates are similar
dyg<-DGEList(counts, group=groups)
dyg<-DGEList(counts, group=grpups)
dyg<-calcNormFactors(dyg, method="TMM")
normd<-cpm(dyg, normalize.lib.sizes=TRUE)
plotMDS(normd, labels=groups) #lets you look for outliers and remove if necessary

##remove genes that have all zero read counts
normd<-normd[ rowSums(normd)!=0, ]

design<-read.csv("YG_design.csv")##load design object for masigpro. row order must be the same as the order of libraries in count matrix
rownames(design)<-design$X #make first column the rownames
design$X<-NULL
d<-make.design.matrix(design, degree=5)

##using the negative binomial options in masigpro, calculate polynomial regressions for each gene
NBp<-p.vector(normd, d, counts=TRUE) #please choose dis. family with care. default for counts is neg. binomial

##TO REMOVE INFLUENTIAL GENES:
NBt<-T.fit(NBp)
influential<-NBt$influ.info
inf.genenames<-colnames(influential)
normd<-normd[!rownames(normd) %in% inf.genenames, ]

##pick k
NBp<-p.vector(normd, d, counts=TRUE)
wss<-(nrow(NBp$SELEC)-1)*sum(apply(NBp$SELEC,2,var))
for (i in 2:15) wss[i]<- sum(kmeans(NBp$SELEC, centers=i, iter.max=20)$withinss)
plot(1:15, wss, type="b")

##estimate m for mfuzz
NBgenes<-ExpressionSet(assayData=NBp$SELEC)
NBgenes.f<-filter.std(NBgenes, min.std=0)
NBgenes.s<-standardise(NBgenes.f)
m1<-mestimate(NBgenes.s)
m1
NBt<-T.fit(NBp)

#load revised see.genes function (separate file, "see.genes.kh")
#Note, I wanted masigpro to cluster all genes regardless if they were differentially expressed under drought, so I modified the underlying code. For most, using the regular see.genes() function will suffice. 
profilesall<-see.genes.kh(NBp$SELEC, edesign=d$edesign, dis=d$dis, cluster.data=1, groups.vector=d$groups.vector, cluster.method="mfuzz", k=9, m=1.06, show.fit=TRUE)

##get D vs W significant genes
sigswd<-get.siggenes(NBt, rsq=0.7, vars="groups", significant.intercept="dummy") #addition of intercept terms includes consideration of control vs. drought intercept values
