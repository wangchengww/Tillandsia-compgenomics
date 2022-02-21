#clean masigpro history
#Karolina Heyduk, University of Georgia, 2016
#modified 2022 (annotations)
#contact: heyduk@hawaii.edu

#load data and libraries
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/")
library(edgeR)
library(maSigPro)
library(mclust)
library(stringr)

trimDat = function (counts) {
  data = read.table(counts,header= T, row.names = 1)
  data = data[,-c(1:5)]
  print(paste0("Number of samples: ", dim(data)[2]))
  print(paste0("Starting number of genes: ", dim(data)[1]))
  data=data[apply(cpm(data),1,function(x){!(mean(x)<1)}),]
  print(paste0("SNumber of genes retaned after trimming: ", dim(data)[1]))
  return(data)
}

counts_trim = trimDat("../counts.Tfas_Tlei_6_timepoints.txt")

#set up edgeR object
groups_list <- data.table::transpose(str_split(colnames(counts_trim), "_"))[c(1,3,4)]
groups <- paste0(groups_list[[1]], "_", groups_list[[2]], "_", groups_list[[3]])
dyg<-DGEList(counts_trim, group=groups)
dyg<-calcNormFactors(dyg, method="TMM")
normd<-cpm(dyg, normalize.lib.sizes=TRUE)
plotMDS(normd, labels=groups) #lets you look for outliers and remove if necessary

##remove genes that have all zero read counts
normd<-normd[ rowSums(normd)!=0, ]

##load design object for masigpro. row order must be the same as the order of libraries in count matrix
design <- data.frame(time=as.factor(c(rep(c(1,2,3,4,5,6)))), 
                     sample=as.factor(c(rep(c(rep(1,6), rep(2,6), rep(3,6), rep(4,6),rep(5,6), rep(6,6)),2))),
                     Tfas = as.factor(c(rep(1, 36), rep(0,36))),
                     Tlei = as.factor(c(rep(0, 36), rep(1,36)))
                     )
rownames(design) <- colnames(counts_trim)
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
sigs <- get.siggenes(NBt, rsq = 0.7, vars = "groups")
suma2Venn(sigs$summary[, c(1:2)])
write.table(sigs$summary$TleivsTfas, file = "Genes_Significant_Tfas-vs-Tlei_0.7pre-trimming.txt", quote = F, sep = "\t", row.names = F)
cluster <- see.genes(sigs$sig.genes$TleivsTfas, edesign = d$edesign, show.fit = T, dis=d$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 9)

save(sigs, d, normd, NBp, NBt, file = "maSigPro_data_run_Tfas-vs-Tlei_0.7-pretrimming.RData")
dat <- load("maSigPro_data_run_Tfas-vs-Tlei_0.7-pretrimming.RData")  

PlotGroups(sigs$sig.genes$TleivsTfas, edesign = , dis = d$dis)
cut <- as.data.frame(cluster$cut)
genes1 <- row.names(cut)[which(cut$`cluster$cut`==1)]
genes1 <- normd[rownames(normd)==genes1, ]
genes1 <- normd[rownames(normd) %in% genes1, ]
PlotGroups(genes1, edesign = design)
##pick k
NBp<-p.vector(normd, d, counts=TRUE)
wss<-(nrow(NBp$SELEC)-1)*sum(apply(NBp$SELEC,2,var))
for (i in 2:15) wss[i]<- sum(kmeans(NBp$SELEC, centers=i, iter.max=20)$withinss)
plot(1:15, wss, type="b")


