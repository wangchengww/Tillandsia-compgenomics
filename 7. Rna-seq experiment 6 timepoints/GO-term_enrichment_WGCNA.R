setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
library("WGCNA")

# Set necessary environment for WGCNA
options(stringsAsFactors = FALSE);
allowWGCNAThreads()

# Load input data
lnames = load(file = "coexpression_input_Tfas_vsd_10c12s.RData");
# Load network data saved in the second part.
lnames = load(file = "coexpression_network_Tfas_vsd_10c12s..RData");
lnames

annot = read.delim("orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt", 
                   header = F, sep = "\t")
annot[,1] <- substring(annot[,1],1,14)
# Match probes in the data set to the probe IDs in the annotation file
probes = colnames(datExpr)
probes2annot = match(probes, annot$V1)
# Get the corresponding Locuis Link IDs
Gene_GOterms = annot$V5[probes2annot];
# $ Choose interesting modules
intModules = c("darkseagreen3", "mediumpurple3", "royalblue")
for (module in intModules){
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modGOs = Gene_GOterms[modGenes];
  # Write them into a file
  fileName = paste("GOterms-", module, ".txt", sep="");
  write.table(as.data.frame(modGOs), file = fileName,
              row.names = FALSE, col.names = FALSE, quote = F)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("GOterms-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10);