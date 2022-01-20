setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
library("WGCNA")
library("stringr")

# Set necessary environment for WGCNA
options(stringsAsFactors = FALSE);
allowWGCNAThreads()

# Load input data
lnames = load(file = "coexpression_input_Tfas_vsd_10c12s.RData");
# Load network data saved in the second part.
lnames = load(file = "coexpression_network_Tfas_vsd_10c12s..RData");
lnames
sign_mods <- read.csv("Table_modules_sign_time_Tfas.csv", header = T, row.names = 1)
modNames <- substring(rownames(sign_mods), 3)

annot = read.delim("orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt", 
                   header = F, sep = "\t")
annot[,1] <- substring(annot[,1],1,14)
# Match probes in the data set to the probe IDs in the annotation file
probes = colnames(datExpr)
probes2annot = match(probes, annot$V1)
# Get the corresponding go terms
Gene_GOterms = annot$V5[probes2annot];
# $ Choose interesting modules
for (module in modNames){
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modGOs = Gene_GOterms[modGenes];
  modGOs <- str_split(modGOs, ",")
  modGOs <- unlist(modGOs)
  # Write them into a file
  fileName = paste("GOterms-", module, ".txt", sep="");
  write.table(as.data.frame(modGOs), file = fileName,
              row.names = FALSE, col.names = FALSE, quote = F)
}
