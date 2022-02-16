# In this script, we make step-by-step coexpression networks following the WGCNA tutorials.
setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
library("WGCNA")
library("ggplot2")
library(reshape2)

# Set necessary environment for WGCNA
options(stringsAsFactors = FALSE);
allowWGCNAThreads()

# Load input data
lnames = load(file = "coexpression_input_Tfas-Tlei_vst.RData");

# Choosing the soft threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Evaluate softpower thresholds for an unsigned network
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softpower_table <- sft[["fitIndices"]]
softpower_table$SFT.R.sq <- format(softpower_table$SFT.R.sq, scientific = F)
write.table(softpower_table, file = "SoftPower_table_Tfas-Tlei_unsigned.txt")
softpower_table <- read.table("SoftPower_table_Tfas-Tlei_unsigned.txt", header = T)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Based on what I have read, one should choose a power that maintains an R^2 as high as
# possible and a mean connectivity between 30 and 100. Therefore, I chose a soft-thresholding
# power of 18, where R^2 is 0.7923 and the mean connectivity is 44
softPower = 7;

# Building the ajacency and Topological Overlap Matrix - this is the co-expression network
adjacency = adjacency(datExpr, power = softPower)

TOM = TOMsimilarity(adjacency);

dissTOM = 1-TOM

# Gene clustering
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
#geneTree2 = hclust(as.dist(dissTOM2), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

table(dynamicMods) 
mean(table(dynamicMods))
write.table(table(dynamicMods), file = "Modules_Tfas-Tlei_soft8_unsigned.txt")

# We end up with 70 modules, the largest one contains 5451 genes, the smallest
# contains 41 genes. Now we display the modules under the dendrogram.
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Merge modules with very similar expression, anything with a correlation higher than 90 %.
# This reduces the number of modules from 70 to 68.
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Threshold for merging modules
MEDissThres = 0.1
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Rename to moduleColors - not really sure what this does
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "coexpression_network_unsigned7_Tfas-Tlei_vst.RData")
lnames = load("coexpression_network_unsigned7_Tfas-Tlei_vst.RData")
lnames = load("coexpression_input_Tfas-Tlei_vst.RData")
lnames
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
module_corr <- cbind(as.data.frame(moduleTraitCor), moduleTraitPvalue)
colnames(module_corr) <- c("corr_species", "corr_sample", "corr_time", 
                           "pvalue_species", "pvalue_sample", "pvalue_time")
write.table(module_corr, file = "Correlation_Modules_to_Traits_Tfas-Tlei_unsigned7_vst.txt", quote = F, sep = "\t")
#Gather some info on significant modules
modules_sign_time <- as.data.frame(moduleTraitCor[moduleTraitPvalue[,1] < 0.05 & moduleTraitPvalue[,3] < 0.05,])
colnames(modules_sign_time) <- c("corr_to_species","corr_to_sample", "corr_to_time")
nr_genes_module <- table(moduleColors)
modules_sign_time$gene_count <- nr_genes_module[substring(rownames(modules_sign_time),3)]
modules_sign_time$pvalue_species <- moduleTraitPvalue[rownames(modules_sign_time),1]
modules_sign_time$pvalue_sample <- moduleTraitPvalue[rownames(modules_sign_time),2]
modules_sign_time$pvalue_time <- moduleTraitPvalue[rownames(modules_sign_time),3]
write.table(modules_sign_time, file = "Modules_sign-Time_Tfas-Tlei_soft7.txt", quote = F, sep = "\t")
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Now I want to focus just on time differences
time = as.data.frame(datTraits$time);
names(time) = "time"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, time, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(time), sep="");
names(GSPvalue) = paste("p.GS.", names(time), sep="")

# Plot for each gene in the module the relationship between module membership 
# (i.e. how strongly the gene belongs to the module) and its significance for time 
# (correlation of expression and time, not actually a p-value)
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for time",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Now create a per-gene summary of the network, containing the module, the correlation of 
# each trait to the time points (significance), the p-value of the significance, and the
# membership and o-value for each module, ranked by their significance for time.
# Create the starting data frame
probes = colnames(datExpr)
geneInfo0 = data.frame(moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for time
modOrder = order(-abs(cor(MEs, time, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                               MMPvalue[, modOrder[mod]]);
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.time));
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "geneInfo_co-expression_Tlei_vsd_10c4s_unsigned8.txt", quote = F, sep = "\t")

# Make gene lists for all modules that are correlated with time so that we can run GOterm enrichment for them
modNames <- substring(rownames(modules_sign_time), 3)
for (module in modNames){
        # Select module probes
        modGenes = rownames(geneInfo[geneInfo$moduleColor == module,])
        # Write them into a file
        fileName = paste("Genes-for-Enrichment_Tfas-Tlei_unsigned7_vst_", module, ".txt", sep="");
        write.table(as.data.frame(modGenes), file = fileName,
                    row.names = FALSE, col.names = FALSE, quote = F)
}
