setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
library('DESeq2')
library("vsn")
library("gplots")
library("RColorBrewer")

# Read in counts and create metadata
counts <- read.delim('counts.Tfas.6_timepoints.txt', header = T, row.names = 1)
metadata <- data.frame(sample=c(rep("Tfas_A",6), rep("Tfas_B",6), rep("Tfas_C",6), 
                                rep("Tfas_D",6),rep("Tfas_E",6), rep("Tfas_F",6)), 
                       time=c(rep(c("0100", "0500", "0900","1300", "1700", "2100"))))

# Normalization with DESeq2 is automatically done with the below function
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~time)
dds <- estimateSizeFactors(dds)
# Filter out lowly expressed genes: every gene where fewer than 8 samples have 5 or more counts is eliminated
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 8
sum(idx)
dds <- dds[idx,]
# Calculate dispersion, DE, etc.
dds <- DESeq(dds)
# Two transformations to reduce dependency between mean and variance
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
# Visualize the relationship of mean and variance of all possible transformations
notAllZero<-(rowSums(counts(dds))>0)
pdf("SdPlots2.pdf", width = 8, height = 10)
meanSdPlot(log2(counts(dds,normalized=T)[notAllZero,]+1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()
