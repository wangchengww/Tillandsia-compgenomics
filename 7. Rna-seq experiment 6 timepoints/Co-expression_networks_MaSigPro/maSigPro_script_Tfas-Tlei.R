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
library("WGCNA")

# For full-gene analysis
counts <- read.table("../counts.Tfas_Tlei_6_timepoints.forR.txt", header = T, row.names = 1)
counts <- counts[,-c(1:5)]

# For exonic-only analysis
counts <- read.table("../counts.Tfas_Tlei_6_timepoints.exons.edited.forR.sum.txt", header = T, row.names = 1)

#set up edgeR object
groups_list <- data.table::transpose(str_split(colnames(counts), "_"))[c(1,3,4)]
groups <- paste0(groups_list[[1]], "_", groups_list[[2]], "_", groups_list[[3]])
dyg<-DGEList(counts, group=groups)
dyg<-calcNormFactors(dyg, method="TMM")
normd <- cpm(dyg, normalized.lib.sizes = T)

# Here we trim lowly-expressed genes. This doesn't change the results much but vastly shortens 
# run time
normd_trim <- normd[rowMeans(normd)>1,]
write.table(normd_trim, file = "counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.txt", sep = "\t", quote = F)

# log transform
normd_log <- cpm(dyg, normalized.lib.sizes = T, log = T)
normd_trim_log <- normd_log[row.names(normd_log) %in% row.names(normd_trim),]
rowmed_Tfas <- apply(normd_trim_log[, c(1:36)],1,median)
rowmed_Tlei <- apply(normd_trim_log[, c(37:72)],1,median)

# Median centering, just for visuals
counts_medcentered_Tfas <- normd_trim_log[, c(1:36)] - rowmed_Tfas
counts_medcentered_Tlei <- normd_trim_log[, c(37:72)] - rowmed_Tlei

normd_trim_log_medcentered <- cbind(counts_medcentered_Tfas, counts_medcentered_Tlei)

write.table(normd_trim_log_medcentered, file = "counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.logtransformed.mediancentered.txt", sep = "\t", quote = F)
write.table(normd_trim_log, file = "counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.logtransformed.txt", sep = "\t", quote = F)

##remove genes that have all zero read counts
normd_trim <- normd_trim[ rowSums(normd_trim)!=0, ]
##load design object for masigpro. row order must be the same as the order of libraries in count matrix
# First experimental group will be the baseline, here it is Tlei
design <- data.frame(time = c(rep(c(1,2,3,4,5,6))), 
                     sample = c(rep(1,6), rep(2,6), rep(3,6), rep(4,6),rep(5,6), rep(6,6), 
                                rep(7,6), rep(8,6), rep(9,6), rep(10,6),rep(11,6), rep(12,6)),
                     Tlei = c(rep(0, 36), rep(1,36)),
                     Tfas = c(rep(1, 36), rep(0,36))
                     )
rownames(design) <- colnames(counts)
d<-make.design.matrix(design, degree=5)

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
sigs2 <- get.siggenes(NBt, rsq = 0.7, vars = "each")
dim(sigs$summary)
genes <- sigs$summary
suma2Venn(sigs$summary[, c(1:2)])
write.table(sigs$summary$TfasvsTlei, file = "Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic.txt", quote = F, sep = "\t", row.names = F)

# Get expressional direction
coefficients <- sigs[["sig.genes"]][["TfasvsTlei"]][["coefficients"]][["betatimexTfas"]]
overexpressedTlei <- coefficients[coefficients > 0] # 268 genes
underexpressedTlei <- coefficients[coefficients < 0] # 236 genes

save(sigs, d, normd, NBp, NBt, file = "maSigPro_data_run_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic.RData")
dat <- load("maSigPro_data_run_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic.RData")  

##pick k
sizeGrWindow(9, 5)
par(mfrow = c(1,1));
cex1 = 0.9;
wss<-(nrow(NBp$SELEC)-1)*sum(apply(NBp$SELEC,2,var))
for (i in 2:15) wss[i]<- sum(kmeans(NBp$SELEC, centers=i, iter.max=20)$withinss)
plot(1:15, wss, type="b")
# The transition is not super clear, it is either 6 or 7. I picked 7.

cluster <- see.genes(sigs$sig.genes$TfasvsTlei, show.fit = T, dis=d$dis,
                     cluster.method="hclust" ,cluster.data = 1, k = 7)
# Cut out genes per cluster
cut <- as.data.frame(cluster$cut)
k = 7
for (i in 1:k){
  print(i)
  genes <- row.names(cut)[which(cut$`cluster$cut`==i)]
  genes <- normd_trim[rownames(normd_trim) %in% genes, ]
  genelist = rownames(genes)
  # Write them into a file
  fileName = paste("Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic-cluster", i, ".txt", sep="");
  write.table(as.data.frame(genelist), file = fileName,
              row.names = FALSE, col.names = FALSE, quote = F)
}
