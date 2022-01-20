# In this script, we evaluate, filter, normalize and trim count data so that it is ready for
# building a coexpression network. We use DESeq2 for normalization, low-expression filtering
# and transformation. We test several filtering thresholds and transformations. Then, using the
# WGCNA package, we set the data up for coexpression analysis.

setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
library('DESeq2')
library("vsn")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("grid")
library("gridExtra")
library("WGCNA")

# Read in counts and create metadata
counts <- read.delim('counts.Tfas.6_timepoints.txt', header = T, row.names = 1)
metadata <- data.frame(sample=c(rep(1,6), rep(2,6), rep(3,6), 
                                rep(4,6),rep(5,6), rep(6,6)), 
                       time=c(rep(c(1,2,3,4,5,6))))

# Turn the counts and metadata into a DESeq object
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~time)
dds <- estimateSizeFactors(dds)
# Filter out lowly expressed genes: I try out several filtering thresholds both on number of counts
# per gene and number of samples. The first command filters out any gene that has less than 10 counts
# in more than 8 samples.
idx_10c8s<- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 8 # 20,681 genes
idx_5c8s <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 8 # 21,964 genes
idx_10c12s <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 12 # 19,983 genes
idx_5c12s <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 12 # 21,106 genes
# Subset the count data
dds_10c8s <- dds[idx_10c8s,]
dds_5c8s <- dds[idx_5c8s,]
dds_10c12s <- dds[idx_10c12s,]
dds_5c12s <- dds[idx_5c12s,]
# Calculate dispersion, DE, etc.
dds_10c8s <- DESeq(dds_10c8s)
dds_5c8s <- DESeq(dds_5c8s)
dds_10c12s <- DESeq(dds_10c12s)
dds_5c12s <- DESeq(dds_5c12s)

# Two transformations to reduce dependency between mean and variance
# 10c8s - vst seems best
rld_10c8s <- rlog(dds_10c8s)
vsd_10c8s <- varianceStabilizingTransformation(dds_10c8s)
# Visualize the relationship of mean and variance of all possible transformations
notAllZero <- (rowSums(counts(dds_10c8s,normalized = T)) > 0 )
msd <- meanSdPlot(log2(counts(dds_10c8s,normalized=T)[notAllZero,]+1)) 
msd2 <- meanSdPlot(assay(rld_10c8s[notAllZero,]))
msd3 <- meanSdPlot(assay(vsd_10c8s[notAllZero,]))
pdf("SdPlots_Tfas_10c8s.pdf", width = 8, height = 10)
msd$gg + scale_y_continuous(limits = c(0, 6))
msd2$gg + scale_y_continuous(limits = c(0, 6))
msd3$gg + scale_y_continuous(limits = c(0, 6))
dev.off()

# 5c8s
rld_5c8s <- rlog(dds_5c8s)
vsd_5c8s <- varianceStabilizingTransformation(dds_5c8s)
# Visualize the relationship of mean and variance of all possible transformations
notAllZero <- (rowSums(counts(dds_5c8s,normalized = T)) > 0 )
msd <- meanSdPlot(log2(counts(dds_5c8s,normalized=T)[notAllZero,]+1)) 
msd2 <- meanSdPlot(assay(rld_5c8s[notAllZero,]))
msd3 <- meanSdPlot(assay(vsd_5c8s[notAllZero,]))
pdf("SdPlots_Tfas_5c8s.pdf", width = 8, height = 10)
msd$gg + scale_y_continuous(limits = c(0, 6))
msd2$gg + scale_y_continuous(limits = c(0, 6))
msd3$gg + scale_y_continuous(limits = c(0, 6))
dev.off()

# 10c12s
rld_10c12s <- rlog(dds_10c12s)
vsd_10c12s <- varianceStabilizingTransformation(dds_10c12s)
# Visualize the relationship of mean and variance of all possible transformations
notAllZero <- (rowSums(counts(dds_10c12s,normalized = T)) > 0 )
msd <- meanSdPlot(log2(counts(dds_10c12s,normalized=T)[notAllZero,]+1)) 
msd2 <- meanSdPlot(assay(rld_10c12s[notAllZero,]))
msd3 <- meanSdPlot(assay(vsd_10c12s[notAllZero,]))
pdf("SdPlots_Tfas_10c12s.pdf", width = 8, height = 10)
msd$gg + scale_y_continuous(limits = c(0, 6))
msd2$gg + scale_y_continuous(limits = c(0, 6))
msd3$gg + scale_y_continuous(limits = c(0, 6))
dev.off()

# 5c12s
rld_5c12s <- rlog(dds_5c12s)
vsd_5c12s <- varianceStabilizingTransformation(dds_5c12s)
# Visualize the relationship of mean and variance of all possible transformations
notAllZero <- (rowSums(counts(dds_5c12s,normalized = T)) > 0 )
msd <- meanSdPlot(log2(counts(dds_5c12s,normalized=T)[notAllZero,]+1)) 
msd2 <- meanSdPlot(assay(rld_5c12s[notAllZero,]))
msd3 <- meanSdPlot(assay(vsd_5c12s[notAllZero,]))
pdf("SdPlots_Tfas_5c12s.pdf", width = 8, height = 10)
msd$gg + scale_y_continuous(limits = c(0, 6))
msd2$gg + scale_y_continuous(limits = c(0, 6))
msd3$gg + scale_y_continuous(limits = c(0, 6))
dev.off()

# Generally, VST seems to be the best transformation regardless of filtering threshold.

# Now, we prepare the data for use in WGCNA
plotPCA(vsd_10c12s, intgroup=c("time", "sample"))
options(stringsAsFactors = FALSE)
# Collect and transpose expression data
datExpr0 <- assay(vsd_10c12s)
datExpr0 <- t(datExpr0)
# test if all are "good genes"
gsg = goodSamplesGenes(assay(vsd_10c12s), verbose = 3);
gsg$allOK
# Visualize outliers
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Sample B and F at 21:00 seem outliers. We remove them manually.
# Plot a line to show the cut
abline(h = 145, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 145, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datTraits <- metadata[-c(12,36),]
rownames(datTraits) <- rownames(datExpr)
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Saving data ready for network construction
save(datExpr, datTraits, file = "coexpression_input_Tfas_vsd_10c12s.RData")
