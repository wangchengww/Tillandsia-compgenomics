#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("BiocManager", "topGO", "stringr")

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

# Load arguments
# 1 is GotoGenes.map, 2 is the subset of genes to test, 3 is name of the output file, 4 is the per-gene file
args <- commandArgs(trailingOnly = TRUE)

 # Run enrichment
geneID2GO <- readMappings(file = args[[1]])
GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)
dup_genes <- read.table(args[[2]])
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
allRes[startsWith(allRes$weight01_pval, "<"),6] <- "1e-30"

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
                          "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE)))


# Add orthogroups
per_gene <- read.delim(args[[4]],
                       sep = "\t", header = F)
orthog.list <- create_orthog_list(go.dataframe, per_gene)
go.dataframe2 = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term,
                           "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE)),
                           "Orthogroups" = as.vector(orthog.list),
                           "Genes" = as.vector(enriched_go_with_my_genes.list)
)
go.dataframe2$adj_pval <- as.numeric(go.dataframe2$adj_pval)

go.dataframe2.significant <- go.dataframe2[go.dataframe2$adj_pval <= 0.05,]

write.table(go.dataframe2.significant, file = args[[3]], sep = "\t",
            quote = F, row.names = F)
