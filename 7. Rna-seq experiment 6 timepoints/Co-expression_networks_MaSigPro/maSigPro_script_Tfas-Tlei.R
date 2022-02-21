#clean masigpro history
#Karolina Heyduk, University of Georgia, 2016
#modified 2022 (annotations)
# Modified by Clara Groot Crego
#contact: clara.groot.crego@univie.ac.at

#load data and libraries
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/")
setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/")
library(edgeR)
library(maSigPro)
library(mclust)
library(stringr)

counts <- read.table("../counts.Tfas_Tlei_6_timepoints.forR.txt", header = T, row.names = 1)
counts <- counts[,-c(1:5)]

#set up edgeR object
groups_list <- data.table::transpose(str_split(colnames(counts_trim), "_"))[c(1,3,4)]
groups <- paste0(groups_list[[1]], "_", groups_list[[2]], "_", groups_list[[3]])
dyg<-DGEList(counts, group=groups)
dyg<-calcNormFactors(dyg, method="TMM")
normd <- cpm(dyg, normalized.lib.sizes = T)
# Here we trim lowly-expressed genes. This doesn't change the results much but vastly shortens 
# run time
normd_trim <- normd[rowMeans(normd)>1,]

##remove genes that have all zero read counts
normd_trim <- normd_trim[ rowSums(normd_trim)!=0, ]
##load design object for masigpro. row order must be the same as the order of libraries in count matrix
design <- data.frame(time = c(rep(c(1,2,3,4,5,6))), 
                     sample = c(rep(c(rep(1,6), rep(2,6), rep(3,6), rep(4,6),rep(5,6), rep(6,6)),2)),
                     Tfas = c(rep(1, 36), rep(0,36)),
                     Tlei = c(rep(0, 36), rep(1,36))
                     )
rownames(design) <- colnames(counts_trim)
d<-make.design.matrix(design, degree=5)

# Needed for curve figures
alt_design <- data.frame(time = c(rep(c(rep(1,6), rep(2,6), rep(3,6), rep(4,6),rep(5,6), rep(6,6)),2)), 
                         sample = c(1:72),
                         Tfas = c(rep(1, 36), rep(0,36)),
                         Tlei = c(rep(0, 36), rep(1,36)))


##using the negative binomial options in masigpro, calculate polynomial regressions for each gene
NBp<-p.vector(normd_trim, d, counts=TRUE) #please choose dis. family with care. default for counts is neg. binomial
NBp$i # returns the number of significant genes

##TO REMOVE INFLUENTIAL GENES: Not a big deal if you pre-trimmed
NBt<-T.fit(NBp)
influential<-NBt$influ.info
inf.genenames<-colnames(influential)
normd<-normd[!rownames(normd) %in% inf.genenames, ]

# Get significant genes 
sigs <- get.siggenes(NBt, rsq = 0.7, vars = "groups")
suma2Venn(sigs$summary[, c(1:2)])
write.table(sigs$summary$TleivsTfas, file = "Genes_Significant_Tfas-vs-Tlei_0.7-trimmed.txt", quote = F, sep = "\t", row.names = F)

save(sigs, d, normd, NBp, NBt, file = "maSigPro_data_run_Tfas-vs-Tlei_0.7-trimmed.RData")
dat <- load("maSigPro_data_run_Tfas-vs-Tlei_0.7-pretrimming.RData")  

cluster <- see.genes(sigs$sig.genes$TleivsTfas, show.fit = T, dis=d$dis,
                     cluster.method="hclust" ,cluster.data = 1, k = 9)
# Cut out genes per cluster
cut <- as.data.frame(cluster$cut)
k = 9
for (i in 1:k){
  print(i)
  genes <- row.names(cut)[which(cut$`cluster$cut`==i)]
  genes <- normd_trim[rownames(normd_trim) %in% genes, ]
  PlotGroups(genes, edesign = alt_design)
  genelist = rownames(genes)
  # Write them into a file
  fileName = paste("Genes_Significant_Tfas-vs-Tlei_0.7-trimmed-cluster", i, ".txt", sep="");
  write.table(as.data.frame(genelist), file = fileName,
              row.names = FALSE, col.names = FALSE, quote = F)
}

pepck = "Tfasc_v1.03128"
pepck_data <- normd[rownames(normd)==pepck, ]
PlotGroups(pepck_data, edesign = alt_design)
##pick k
NBp<-p.vector(normd, d, counts=TRUE)
wss<-(nrow(NBp$SELEC)-1)*sum(apply(NBp$SELEC,2,var))
for (i in 2:15) wss[i]<- sum(kmeans(NBp$SELEC, centers=i, iter.max=20)$withinss)
plot(1:15, wss, type="b")


