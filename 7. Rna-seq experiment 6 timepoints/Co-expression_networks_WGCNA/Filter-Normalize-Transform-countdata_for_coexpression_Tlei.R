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
counts <- read.delim('counts.Tlei.6_timepoints.txt', header = T, row.names = 1)
metadata <- data.frame(sample=c(rep(1,6), rep(2,6), rep(3,6), 
                                rep(4,6),rep(5,6), rep(6,6)), 
                       time=c(rep(c(1,2,3,4,5,6))))
metadata$sample <- as.factor(metadata$sample)
metadata$time <- as.factor(metadata$time)

# Turn the counts and metadata into a DESeq object
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~time)
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
idx_10c4s<- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 4 # 21,559 genes
sum(idx_10c4s)
idx_5c4s <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 4 # 22,808 genes
sum(idx_5c4s)
idx_10c6s <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 6 # 21,056 genes
sum(idx_10c6s)
idx_5c6s <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 6 # 22,244 genes
sum(idx_5c6s)
# Subset the count data
dds_10c4s <- dds[idx_10c4s,]

write.table(counts(dds_10c4s, normalized=TRUE), file = "counts.Tlei.6_timepoints.filtr-normalized_DESEQ2.txt")

# The largest difference between filtering sets is 1,945 genes, which is only 6 % of the total
# amount of genes. I decided to carry on with the stringent count criterium but non-stringent sample
# criterium (less than 10 counts in more than 32 samples, which includes 21,818 genes)
# Calculate dispersion, DE, etc.
dds_10c4s <- DESeq(dds_10c4s)

# Two transformations to reduce dependency between mean and variance
# 10c8s - vst seems best
rld_10c4s <- rlog(dds_10c4s)
vsd_10c4s <- varianceStabilizingTransformation(dds_10c4s)
# Visualize the relationship of mean and variance of all possible transformations
´
# Generally, VST seems to be the best transformation regardless of filtering threshold.

# Now, we prepare the data for use in WGCNA
plotPCA(vsd_10c4s, intgroup=c("time", "sample"))
options(stringsAsFactors = FALSE)
# Collect and transpose expression data
datExpr0 <- assay(vsd_10c4s)
write.table(datExpr0, file = "counts.Tfas.6_timepoints.filtr-transform-normalized_DESEQ2.txt")
datExpr0 <- t(datExpr0)
# test if all are "good genes"
gsg = goodSamplesGenes(assay(vsd_10c4s), verbose = 3);
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
abline(h = 160, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datTraits <- metadata
rownames(datTraits) <- rownames(datExpr)
datTraits$sample <- as.numeric(datTraits$sample)
datTraits$time <- as.numeric(datTraits$time)
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Saving data ready for network construction
save(datExpr, datTraits, file = "coexpression_input_Tlei_vsd_10c4s.RData")
