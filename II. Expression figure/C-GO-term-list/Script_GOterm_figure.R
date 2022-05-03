setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/C-GO-term-list/')
library("GOplot")
library(stringr)
go <- read.delim("GO-term_enrichment_mod_TLEI-REF_exonic.txt", header = T, sep = "\t")
ortho_info <- read.delim("orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.blastandsearch_noAcom.txt",
                          sep = "\t")
colnames(ortho_info) <- c("gene_ID", "chr", "start", "end", "GOterm", "Description", "orthogroup", "count_Acom", 
                          "count_Tfas", "count_Tlei")
genelist <- unlist(str_split(go$Genes, ", "))
genes <- subset(ortho_info, gene_ID %in% genelist)

genes <- genes[,c(1,9,10)]
genes$difference <- genes$count_Tfas - genes$count_Tlei
genes <- genes[,c(1,4)]
colnames(genes) <- c("ID", "logFC")
circ <- circle_dat(go, genes)
GOBar(subset(circ, category == 'BP'))
