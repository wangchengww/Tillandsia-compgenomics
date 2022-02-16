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
counts <- read.delim('counts.Tfas_Tlei_6_timepoints.txt', header = T, row.names = 1)
counts <- counts[,-c(1:5)]
metadata <- data.frame(species=c(rep(0,36), rep(1,36)),
                       sample=c(rep(1,6), rep(2,6), rep(3,6), 
                                rep(4,6),rep(5,6), rep(6,6)), 
                       time=c(rep(c(1,2,3,4,5,6))))
metadata$species <- as.factor(metadata$species)
metadata$sample <- as.factor(metadata$sample)
metadata$time <- as.factor(metadata$time)

# Turn the counts and metadata into a DESeq object
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~species*time)
dds <- estimateSizeFactors(dds)
# save normalized counts
counts.norm <- counts(dds, normalized=TRUE)
# Filter out lowly expressed genes: I try out several filtering thresholds both on number of counts
# per gene and number of samples. The first command keeps any gene that has more than 10 counts
# in more than 4 samples, meaning that every gene with more than 88 % of samples with an expression
# under 10 will be filtered out. In the second command, I reduce the count threshold to 5, which is 
# less stringent. In the third command, I keep the count threshold the same (10), but this
# time any gene with more than 83 % of of samples with an expression under 10 will be removed, 
# so it is again more stringent.
idx<- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 7 # 22,381 genes
sum(idx)
# Subset the count data
dds <- dds[idx,]

write.table(counts(dds, normalized=TRUE), file = "counts.Tfas_Tlei_6_timepoints.filtr-normalized_DESEQ2.txt")
t <- read.delim("counts.Tfas_Tlei_6_timepoints.filtr-normalized_DESEQ2.txt", sep = " ")
tfas <- t[, c(1:36)]
write.table(tfas, file = "counts.Tfas.6_timepoints.filtr-JOINTLY-normalized_DESEQ2.txt", sep = "\t", quote = F)
tlei <- t[, c(37:72)]
write.table(tlei, file = "counts.Tlei.6_timepoints.filtr-JOINTLY-normalized_DESEQ2.txt", sep = "\t", quote = F)

dds <- DESeq(dds)

# Two transformations to reduce dependency between mean and variance
# 10c8s - vst seems best
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
# Visualize the relationship of mean and variance of all possible transformations
notAllZero <- (rowSums(counts(dds,normalized = T)) > 0 )
msd <- meanSdPlot(log2(counts(dds,normalized=T)[notAllZero,]+1)) 
msd2 <- meanSdPlot(assay(rld[notAllZero,]))
msd3 <- meanSdPlot(assay(vsd[notAllZero,]))
pdf("SdPlots_Tfas-Tlei.pdf", width = 8, height = 10)
msd$gg + scale_y_continuous(limits = c(0, 6))
msd2$gg + scale_y_continuous(limits = c(0, 6))
msd3$gg + scale_y_continuous(limits = c(0, 6))
dev.off()
# Generally, VST seems to be the best transformation regardless of filtering threshold.

# Now, we prepare the data for use in WGCNA
options(stringsAsFactors = FALSE)
# Collect and transpose expression data
datExpr0 <- assay(vsd)
datExpr0 <- t(datExpr0)
# test if all are "good genes"
gsg = goodSamplesGenes(assay(vsd), verbose = 3);
gsg$allOK
# Visualize outliers
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Sample Tfas_B and Tfas_F at 21:00 seem outliers. We remove them manually.
# Plot a line to show the cut
abline(h = 170, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 170, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust!=0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datTraits <- metadata[-c(12,36),]
rownames(datTraits) <- rownames(datExpr)
datTraits$sample <- as.numeric(datTraits$sample)
datTraits$time <- as.numeric(datTraits$time)
datTraits$species <- as.numeric(datTraits$species)
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Saving data ready for network construction
save(datExpr, datTraits, file = "coexpression_input_Tfas-Tlei_vst.RData")
