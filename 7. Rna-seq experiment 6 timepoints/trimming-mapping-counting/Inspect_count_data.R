setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/")
library("edgeR")
library("RUVSeq")
library("HTSFilter")
library("devtools")
library("feathers")
library("RColorBrewer")

#Function to trim the data
trimDat = function (counts) {
  data = read.table(counts,header= T, row.names = 1)
  #data = data[,c(6:77)]
  print(paste0("Number of samples: ", dim(data)[2]))
  print(paste0("Starting number of genes: ", dim(data)[1]))
  data=data[apply(cpm(data),1,function(x){!(mean(x)<1)}),]
  print(paste0("SNumber of genes retaned after trimming: ", dim(data)[1]))
  return(data)
}

data <- read.table("counts.Tfas_Tlei_6_timepoints.txt",header= T, row.names = 1) 
data <- data[,6:77]

# Look at count distribution per sample
par(mfrow=c(2,3))
for (col in 1:ncol(data)) {
  hist(data[,col][sample(1:11900,100)], 
       breaks = 50,
       xlab = "read count",
       main = paste0("Sample number ",col))
}

# Look at log(cpm) distribution across genes? and expression variance across log(cpm)
par(mfrow=c(1,2))
hist(log(cpm(data)),
     breaks = 1000,
     main = "Frequency of log(cpm)",
     xlab = "Expression (log(cpm))")
plot(apply(log(cpm(data)), 1, mean), apply(log(cpm(data)), 1, var), 
     xlab="Mean expression (log(cpm))",
     ylab="Variance of expression log(cpm)", 
     main = "Expression Variance")

# Trim low-expression genes and remake plots
 data_trim = trimDat("counts.Tfas_Tlei_6_timepoints.txt")
 write.table(data_trim, file ="counts.Tfas_Tlei_6_timepoints.filtered.txt", 
             sep = "\t", quote = F)
data_trim <- read.delim("counts.Tfas_Tlei_6_timepoints.filtered.txt")

hist(log(cpm(data_trim)),
     breaks = 1000,
     main = "Frequency of log(cpm)",
     xlab = "Expression (log(cpm))")
plot(apply(log(cpm(data_trim)), 1, mean),apply(log(cpm(data_trim)), 1, var), 
     xlab="Mean expression (log(cpm))",
     ylab="Variance of expression log(cpm)", 
     main = "Expression Variance")

# Plot PCA
groupGLM <- factor(c(rep("Tfas.A",6),rep("Tfas.B",6),rep("Tfas.C",6), # Populations
                     rep("Tfas.D",6),rep("Tfas.E",6),rep("Tfas.F",6),
                     rep("Tlei.A",6),rep("Tlei.C",6),rep("Tlei.D",6),
                     rep("Tlei.E",6),rep("Tlei.F",6),rep("Tlei.G",6)))

y <- DGEList(counts=data_trim,group=groupGLM)
samples <- c(rep("#E41A1C",6),rep("#377EB8",6), rep("#4DAF4A",6), #define colors for ecotype pairs
          rep("#984EA3",6),rep("#FF7F00",6), rep("#FFFF33",6),
          rep("#E41A1C",6),rep("#4DAF4A",6),rep("#984EA3",6),
          rep("#FF7F00",6),rep("#FFFF33",6),rep("#A65628",6))
times <- c(rep(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"),12))
species <- c(rep(17,36),rep(19,36)) #define symbols for ecotypes
#Plot PCA, change k for visualizing different numbers of components
par(mfrow=c(1,1))
plotPCA(y$counts,col=samples,  k=5, cex = 1.7, pch=species, main = "PCA of transcripts counts", labels = F)
plotPCA(y$counts,col=times,  k=2, cex = 1.7, pch=species, main = "PCA of transcripts counts", labels = F)

#--------------------------- Inspection of per-species dataset - T.fasciculata ----------------------

data_tfas <- data[,c(1:36)]
write.table(data_tfas, "counts.Tfas.6_timepoints.txt", sep = "\t", quote = F)
# Trim low-expression genes and remake plots
# data_tfas_trim = trimDat("counts.Tfas.6_timepoints.txt")
# write.table(data_tfas_trim, file ="counts.Tfas.6_timepoints.filtered.txt", 
#              sep = "\t", quote = F)
data_tfas_trim <- read.delim("counts.Tfas.6_timepoints.filtered.txt")

# Plot PCA
groupGLM <- factor(c(rep("Tfas.A",6),rep("Tfas.B",6),rep("Tfas.C",6), # Populations
                     rep("Tfas.D",6),rep("Tfas.E",6),rep("Tfas.F",6)))

y <- DGEList(counts=data_tfas_trim,group=groupGLM)
samples <- c(rep("#E41A1C",6),rep("#377EB8",6), rep("#4DAF4A",6), #define colors for ecotype pairs
             rep("#984EA3",6),rep("#FF7F00",6), rep("#FFFF33",6))
times <- c(rep(c(15,16,17,18,19,8),6)) #define symbols for ecotypes
#Plot PCA, change k for visualizing different numbers of components
par(mfrow=c(1,1))
plotPCA(y$counts,col=samples,  k=6, cex = 1.7, pch=times, main = "PCA of transcripts counts - T. fasciculata", labels = F)

#----------------------- Inspection of per-species dataset - T.leiboldiana-------------------------

data_tlei <- data[,c(37:72)]
write.table(data_tlei, "counts.Tlei.6_timepoints.txt", sep = "\t", quote = F)
# Trim low-expression genes and remake plots
# data_tlei_trim = trimDat("counts.Tlei.6_timepoints.txt")
# write.table(data_tlei_trim, file ="counts.Tlei.6_timepoints.filtered.txt", 
#              sep = "\t", quote = F)
data_tlei_trim <- read.delim("counts.Tlei.6_timepoints.filtered.txt")

# Plot PCA
groupGLM <- factor(c(rep("Tlei.A",6),rep("Tlei.C",6),rep("Tlei.D",6), # Populations
                     rep("Tlei.E",6),rep("Tlei.F",6),rep("Tlei.G",6)))

y <- DGEList(counts=data_tlei_trim,group=groupGLM)
samples <- c(rep("#E41A1C",6),rep("#4DAF4A",6),rep("#984EA3",6),
             rep("#FF7F00",6),rep("#FFFF33",6),rep("#A65628",6))
times <- c(rep(c(15,16,17,18,19,8),6)) #define symbols for ecotypes
#Plot PCA, change k for visualizing different numbers of components
par(mfrow=c(1,1))
plotPCA(y$counts,col=samples,  k=4, cex = 1.7, pch=times, main = "PCA of transcripts counts - T. leiboldiana", labels = F)


#------------------------- Inspection of count set of ORTHOLOGOUS GENES only ------------------------

# data_trim = trimDat("counts.Tfas_Tlei_6_timepoints.orthologs.txt")
# write.table(data_trim, file ="counts.Tfas_Tlei_6_timepoints.orthologs.filtered.txt", 
#             sep = "\t", quote = F)
data <- read.delim("counts.Tfas_Tlei_6_timepoints.orthologs.txt", header = T, sep = '\t', row.names=1)
data_trim <- read.delim("counts.Tfas_Tlei_6_timepoints.orthologs.filtered.txt", header = T, sep = "\t")

hist(log(cpm(data_trim)),
     breaks = 1000,
     main = "Frequency of log(cpm)",
     xlab = "Expression (log(cpm))")
plot(apply(log(cpm(data_trim)), 1, mean),apply(log(cpm(data_trim)), 1, var), 
     xlab="Mean expression (log(cpm))",
     ylab="Variance of expression log(cpm)", 
     main = "Expression Variance")

# Plot PCA
groupGLM <- factor(c(rep("Tfas.A",6),rep("Tfas.B",6),rep("Tfas.C",6), # Populations
                     rep("Tfas.D",6),rep("Tfas.E",6),rep("Tfas.F",6),
                     rep("Tlei.A",6),rep("Tlei.C",6),rep("Tlei.D",6),
                     rep("Tlei.E",6),rep("Tlei.F",6),rep("Tlei.G",6)))

y <- DGEList(counts=data_trim,group=groupGLM)
samples <- c(rep("#E41A1C",6),rep("#377EB8",6), rep("#4DAF4A",6), #define colors for ecotype pairs
             rep("#984EA3",6),rep("#FF7F00",6), rep("#FFFF33",6),
             rep("#E41A1C",6),rep("#4DAF4A",6),rep("#984EA3",6),
             rep("#FF7F00",6),rep("#FFFF33",6),rep("#A65628",6))
times <- c(rep(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"),12))
species <- c(rep(17,36),rep(19,36)) #define symbols for ecotypes
#Plot PCA, change k for visualizing different numbers of components
par(mfrow=c(1,1))
plotPCA(y$counts,col=samples,  k=5, cex = 1.7, pch=species, main = "PCA of transcripts counts", labels = F)
plotPCA(y$counts,col=times,  k=2, cex = 1.7, pch=species, main = "PCA of transcripts counts", labels = F)


# Orthologous genes - T. fasciculata samples only

data_tfas <- data[,c(6:41)]
write.table(data_tfas, "counts.Tfas.6_timepoints.orthologs.txt", sep = "\t", quote = F)
# Trim low-expression genes and remake plots
data_tfas_trim = trimDat("counts.Tfas.6_timepoints.orthologs.txt")
write.table(data_tfas_trim, file ="counts.Tfas.6_timepoints.orthologs.filtered.txt", 
              sep = "\t", quote = F)
data_tfas_trim <- read.delim("counts.Tfas.6_timepoints.orthologs.filtered.txt")
# Plot PCA
groupGLM <- factor(c(rep("Tfas.A",6),rep("Tfas.B",6),rep("Tfas.C",6), # Populations
                     rep("Tfas.D",6),rep("Tfas.E",6),rep("Tfas.F",6)))

y <- DGEList(counts=data_tfas_trim,group=groupGLM)
samples <- c(rep("#E41A1C",6),rep("#377EB8",6), rep("#4DAF4A",6), #define colors for ecotype pairs
             rep("#984EA3",6),rep("#FF7F00",6), rep("#FFFF33",6))
times <- c(rep(c(15,16,17,18,19,8),6)) #define symbols for ecotypes
#Plot PCA, change k for visualizing different numbers of components
par(mfrow=c(1,1))
plotPCA(y$counts,col=samples,  k=4, cex = 1.7, pch=times, main = "PCA of transcripts counts - T. fasciculata", labels = F)

# Orthologous genes - T. leiboldiana samples only

data_tlei <- data[,c(42:77)]
write.table(data_tlei, "counts.Tlei.6_timepoints.orthologs.txt", sep = "\t", quote = F)
# Trim low-expression genes and remake plots
data_tlei_trim = trimDat("counts.Tlei.6_timepoints.orthologs.txt")
write.table(data_tlei_trim, file ="counts.Tlei.6_timepoints.orthologs.filtered.txt", 
            sep = "\t", quote = F)
data_tlei_trim <- read.delim("counts.Tlei.6_timepoints.orthologs.filtered.txt")

# Plot PCA
groupGLM <- factor(c(rep("Tlei.A",6),rep("Tlei.C",6),rep("Tlei.D",6), # Populations
                     rep("Tlei.E",6),rep("Tlei.F",6),rep("Tlei.G",6)))

y <- DGEList(counts=data_tlei_trim,group=groupGLM)
samples <- c(rep("#E41A1C",6),rep("#4DAF4A",6),rep("#984EA3",6),
             rep("#FF7F00",6),rep("#FFFF33",6),rep("#A65628",6))
times <- c(rep(c(15,16,17,18,19,8),6)) #define symbols for ecotypes
#Plot PCA, change k for visualizing different numbers of components
par(mfrow=c(1,1))
plotPCA(y$counts,col=samples,  k=4, cex = 1.7, pch=times, main = "PCA of transcripts counts - T. leiboldiana", labels = F)

