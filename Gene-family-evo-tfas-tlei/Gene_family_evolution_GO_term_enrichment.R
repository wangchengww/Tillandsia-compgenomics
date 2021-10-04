# Studying differences in gene family sizes between Tfas and Tlei

library(ggplot2)
library(dplyr)

setwd("Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
counts <- read.table("orthogroups_Tfas_Tlei_Acom.counts.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.blastandsearch.txt", sep = '\t')
colnames(counts) <- c("og_id", "Acom", "Tfas", "Tlei")
per_gene <- read.delim("orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.blastandsearch_noAcom.txt", 
                       sep = "\t", header = F)

# Filter out unique orthogroups
counts_Tfas_Tlei <- counts[counts$Tfas != 0 & counts$Tlei != 0,]
ggplot(counts_Tfas_Tlei, aes(x=Tfas, y=Tlei)) + geom_point(size = .5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Per-species gene counts in multi-copy orthogroups") +
  ylab(label = "T. leiboldiana") +
  xlab(label = "T. fasciculata") +
  theme_bw()

# Add colour gradient to show number of datapoints with the same count combination
counts_Tfas_Tlei_multi <- counts_Tfas_Tlei[!(counts_Tfas_Tlei$Tfas == 1 & counts_Tfas_Tlei$Tlei == 1),]
ggplot(counts_Tfas_Tlei_multi) + geom_hex(aes(Tfas, Tlei, fill = stat(log(count))), bins = 100) +
  labs(title = "Per-species gene counts in multi-copy orthogroups") +
  ylab(label = "T. leiboldiana") +
  xlab(label = "T. fasciculata") +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

counts_more_Tlei <- counts_Tfas_Tlei_multi[(counts_Tfas_Tlei_multi$Tfas < counts_Tfas_Tlei_multi$Tlei),]
counts_more_Tfas <- counts_Tfas_Tlei_multi[(counts_Tfas_Tlei_multi$Tfas > counts_Tfas_Tlei_multi$Tlei),]

write.table(counts_Tfas_Tlei_multi, file = "orthogroup_selection_multicopy_for_GO_term_all.txt", sep = "\t", quote = F, row.names = F)
write.table(counts_more_Tfas, file = "orthogroup_selection_multicopy_larger_in_Tfas.txt", sep = "\t", quote = F, row.names = F)
write.table(counts_more_Tlei, file = "orthogroup_selection_multicopy_larger_in_Tlei.txt", sep = "\t", quote = F, row.names = F)

### Log-ratio test
counts_Tfas_Tlei$logratio <- log(counts_Tfas_Tlei$Tfas/counts_Tfas_Tlei$Tlei)
mean_logratio <- mean(counts_Tfas_Tlei$logratio) # 0.0203
counts_Tfas_Tlei$corr_logratio <- counts_Tfas_Tlei$logratio - mean_logratio


top2percent_Tfas_larger <- counts_Tfas_Tlei %>%
   arrange(desc(corr_logratio)) %>%
   filter(corr_logratio >= quantile(corr_logratio, .98))
top2percent_Tlei_larger <- counts_Tfas_Tlei %>%
   arrange((corr_logratio)) %>%
   filter(corr_logratio <= quantile(corr_logratio, .02))

top2percent_Tfas_larger_multi <- counts_Tfas_Tlei_multi %>%
   arrange(desc(corr_logratio)) %>% 
   filter(corr_logratio >= quantile(corr_logratio, .98))
top2percent_Tlei_larger_multi <- counts_Tfas_Tlei_multi %>%
   arrange((corr_logratio)) %>% 
   filter(corr_logratio <= quantile(corr_logratio, .02))

write.table(top2percent_Tfas_larger_multi, file = "Top2percent_multicopy_orthogroups_LogRatio_larger_in_Tfas.txt", sep = "\t", quote = F, row.names = F)
write.table(top2percent_Tlei_larger_multi, file = "Top2percent_multicopy_orthogroups_LogRatio_larger_in_Tlei.txt", sep = "\t", quote = F, row.names = F)


# Because of the very high occurrence of 1:1 orthogroups, anything deviating from 
# this will already be in the top 2 % (includes all duplications, also 2:1). This shows that any
# form of copy variance is in fact already "deviating" as far as we can tell.

### GO term enrichment on all multi-copy genes ###

# Installation
# if (!requireNamespace("BiocManager", quietly=TRUE)) + install.packages("BiocManager")
# BiocManager::install()
# if (!requireNamespace("BiocManager", quietly=TRUE)) + install.packages("BiocManager")
# BiocManager::install("topGO")
# if (!requireNamespace("BiocManager", quietly=TRUE)) + install.packages("BiocManager")
#  BiocManager::install("ALL")
 
# Load packages
library(topGO)
library(stringr)

# Load functions
 change_names <- function(data, name_list){
   colnames(data) <- name_list
   return(data)
 }
 
 rename <- function(table, geneNames){
   names(table) <- geneNames
   return(table)
 }
 
 attach_enriched_go_genes <- function(enriched_go_with_my_genes){
   enriched_go_with_my_genes.list = c()
   for (i in 1:length(enriched_go_with_my_genes)){
     enriched_go_with_my_genes.list = c(enriched_go_with_my_genes.list, enriched_go_with_my_genes[[i]])
   }
   return(enriched_go_with_my_genes.list)
 }
 
 create_orthog_list <- function(go.dataframe, per_gene){
   genes <- str_split(go.dataframe[,4], ", ")
   orthog.l <- vector(mode = "list", length = 300)
   l = 1
   for (i in genes){
     orthog <- vector()
     for (j in i){
       new <- as.character(per_gene[per_gene$V1==j,7])
       orthog <- c(orthog, new)
     }
     orthog_uniq <- unique(orthog)
     orthog_string <- paste(orthog_uniq, collapse = ', ')
     orthog.l[[l]] <- orthog_string
     l = l+1
     orthog <- c()
   }
   names(orthog.l) <- c(go.dataframe$ID)
   orthog.list = attach_enriched_go_genes(orthog.l)
   return(orthog.list)
 }

 # Run enrichment for all multicopy genes
geneID2GO <- readMappings(file = "genes_to_GO.map")
GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)
dup_genes <- read.table("dup_genes_Tfas-Tlei_ID.txt")
geneList <- factor(as.integer(geneNames %in% dup_genes$V1))
names(geneList) <- geneNames
str(geneList)

name_list = c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
table = as.factor(geneNames) %in% dup_genes$V1
int_table = as.integer(table)
int_fac_table = factor(int_table)
fac_table = rename(table = int_fac_table, geneNames = geneNames)

GOdata.BP = new("topGOdata", ontology = "BP", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.MF = new("topGOdata", ontology = "MF", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.CC = new("topGOdata", ontology = "CC", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultWeight01.BP = runTest(GOdata.BP, statistic = "fisher")
resultWeight01.MF = runTest(GOdata.MF, statistic = "fisher")
resultWeight01.CC = runTest(GOdata.CC, statistic = "fisher")

allRes.BP1 = GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.BP2 = cbind(allRes.BP1,"BP")
allRes.BP = change_names(data = allRes.BP2, name_list = name_list)

allRes.MF1 = GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.MF2 = cbind(allRes.MF1,"MF")
allRes.MF = change_names(data = allRes.MF2, name_list = name_list)

allRes.CC1 = GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.CC2 = cbind(allRes.CC1,"CC")
allRes.CC = change_names(data = allRes.CC2, name_list = name_list)

allRes1 = rbind(allRes.BP,allRes.MF)
allRes = rbind(allRes1, allRes.CC)

allGO.BP = genesInTerm(GOdata.BP)
allGO.MF = genesInTerm(GOdata.MF)
allGO.CC = genesInTerm(GOdata.CC)
allGO = c(allGO.BP, allGO.MF, allGO.CC)

# Create final table
SAM_ANOTATION = lapply(allGO,function(x) x[x %in%  dup_genes$V1])
enriched_go_with_my_genes = lapply(SAM_ANOTATION[allRes[,1]], paste0, collapse = ", ")
enriched_go_with_my_genes.list = attach_enriched_go_genes(enriched_go_with_my_genes)
go.dataframe = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term, 
                          "Genes" = as.vector(enriched_go_with_my_genes.list), 
                          "adj_pval" = (sub(",", ".", allRes$weight01_pval, fixed = TRUE)))
go.dataframe[101,5] <- "1e-30"
go.dataframe$adj_pval <- as.numeric(go.dataframe$adj_pval)

# Add orthogroups
orthog.list <- create_orthog_list(go.dataframe, per_gene)
go.dataframe2 = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term, 
                           "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE)),
                           "Orthogroups" = as.vector(orthog.list),
                           "Genes" = as.vector(enriched_go_with_my_genes.list) 
)
new_DF <- go.dataframe2[rowSums(is.na(go.dataframe2)) > 0,]
go.dataframe2[101,5] <- "1e-30"
go.dataframe2$adj_pval <- as.numeric(go.dataframe2$adj_pval)

write.table(go.dataframe2, file = "Enriched_GO_terms_Tfas-Tlei_ALL_multicopy_genes.txt", sep = "\t", 
            quote = F, row.names = F)

### GO term enrichment on orthogroups larger in Tfas ###

dup_genes <- read.table("dup_genes_larger_in_Tfas.txt")
geneList <- factor(as.integer(geneNames %in% dup_genes$V1))
names(geneList) <- geneNames
str(geneList)

name_list = c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
table = as.factor(geneNames) %in% dup_genes$V1
int_table = as.integer(table)
int_fac_table = factor(int_table)
fac_table = rename(table = int_fac_table, geneNames = geneNames)

GOdata.BP = new("topGOdata", ontology = "BP", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.MF = new("topGOdata", ontology = "MF", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.CC = new("topGOdata", ontology = "CC", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultWeight01.BP = runTest(GOdata.BP, statistic = "fisher")
resultWeight01.MF = runTest(GOdata.MF, statistic = "fisher")
resultWeight01.CC = runTest(GOdata.CC, statistic = "fisher")

allRes.BP1 = GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.BP2 = cbind(allRes.BP1,"BP")
allRes.BP = change_names(data = allRes.BP2, name_list = name_list)

allRes.MF1 = GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.MF2 = cbind(allRes.MF1,"MF")
allRes.MF = change_names(data = allRes.MF2, name_list = name_list)

allRes.CC1 = GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.CC2 = cbind(allRes.CC1,"CC")
allRes.CC = change_names(data = allRes.CC2, name_list = name_list)

allRes1 = rbind(allRes.BP,allRes.MF)
allRes = rbind(allRes1, allRes.CC)

allGO.BP = genesInTerm(GOdata.BP)
allGO.MF = genesInTerm(GOdata.MF)
allGO.CC = genesInTerm(GOdata.CC)
allGO = c(allGO.BP, allGO.MF, allGO.CC)

SAM_ANOTATION = lapply(allGO,function(x) x[x %in%  dup_genes$V1])
enriched_go_with_my_genes = lapply(SAM_ANOTATION[allRes[,1]], paste0, collapse = ", ")
enriched_go_with_my_genes.list = attach_enriched_go_genes(enriched_go_with_my_genes)
go.dataframe = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term, 
                          "Genes" = as.vector(enriched_go_with_my_genes.list), 
                          "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE)))
go.dataframe[101,5] <- "1e-30"
go.dataframe$adj_pval <- as.numeric(go.dataframe$adj_pval)

# Add orthogroups
orthog.list <- create_orthog_list(go.dataframe, per_gene)
go.dataframe2 = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term, 
                           "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE)),
                           "Orthogroups" = as.vector(orthog.list),
                           "Genes" = as.vector(enriched_go_with_my_genes.list) 
)
new_DF <- go.dataframe2[rowSums(is.na(go.dataframe2)) > 0,]
go.dataframe2[101,5] <- "1e-30"
go.dataframe2$adj_pval <- as.numeric(go.dataframe2$adj_pval)

write.table(go.dataframe2, file = "Enriched_GO_terms_multicopy_genes_larger-in-Tfas.txt", sep = "\t", 
            quote = F, row.names = F)

### GO term enrichment on orthogroups larger in Tlei ###

dup_genes <- read.table("dup_genes_larger_in_Tlei.txt")
geneList <- factor(as.integer(geneNames %in% dup_genes$V1))
names(geneList) <- geneNames
str(geneList)

name_list = c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
table = as.factor(geneNames) %in% dup_genes$V1
int_table = as.integer(table)
int_fac_table = factor(int_table)
fac_table = rename(table = int_fac_table, geneNames = geneNames)

GOdata.BP = new("topGOdata", ontology = "BP", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.MF = new("topGOdata", ontology = "MF", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.CC = new("topGOdata", ontology = "CC", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultWeight01.BP = runTest(GOdata.BP, statistic = "fisher")
resultWeight01.MF = runTest(GOdata.MF, statistic = "fisher")
resultWeight01.CC = runTest(GOdata.CC, statistic = "fisher")

allRes.BP1 = GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.BP2 = cbind(allRes.BP1,"BP")
allRes.BP = change_names(data = allRes.BP2, name_list = name_list)

allRes.MF1 = GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.MF2 = cbind(allRes.MF1,"MF")
allRes.MF = change_names(data = allRes.MF2, name_list = name_list)

allRes.CC1 = GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 100, numChar=1000)
allRes.CC2 = cbind(allRes.CC1,"CC")
allRes.CC = change_names(data = allRes.CC2, name_list = name_list)

allRes1 = rbind(allRes.BP,allRes.MF)
allRes = rbind(allRes1, allRes.CC)

allGO.BP = genesInTerm(GOdata.BP)
allGO.MF = genesInTerm(GOdata.MF)
allGO.CC = genesInTerm(GOdata.CC)
allGO = c(allGO.BP, allGO.MF, allGO.CC)

SAM_ANOTATION = lapply(allGO,function(x) x[x %in%  dup_genes$V1])
enriched_go_with_my_genes = lapply(SAM_ANOTATION[allRes[,1]], paste0, collapse = ", ")
enriched_go_with_my_genes.list = attach_enriched_go_genes(enriched_go_with_my_genes)

go.dataframe = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term, 
                          "Genes" = as.vector(enriched_go_with_my_genes.list), 
                          "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE)))

orthog.list <- create_orthog_list(go.dataframe, per_gene)

go.dataframe2 = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term, 
                           "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE)),
                         "Orthogroups" = as.vector(orthog.list),
                         "Genes" = as.vector(enriched_go_with_my_genes.list) 
                         )
new_DF <- go.dataframe2[rowSums(is.na(go.dataframe2)) > 0,]
go.dataframe2$adj_pval <- as.numeric(go.dataframe2$adj_pval)

write.table(go.dataframe2, file = "Enriched_GO_terms_multicopy_genes_larger-in-Tlei.txt", sep = "\t", 
            quote = F, row.names = F)
