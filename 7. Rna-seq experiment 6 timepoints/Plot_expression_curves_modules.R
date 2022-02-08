setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/GO-term_Enrichment_Modules/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/GO-term_Enrichment_Modules/')
library("ggplot2")
library("reshape2")
library("stringr")

counts <- read.table("../counts.Tlei.6_timepoints.filtr-normalized_DESEQ2.logtransformed.txt", header = T, row.names = 1)
#counts_log <- log(counts+1)
#write.table(counts_log, file = "../counts.Tfas.6_timepoints.filtr-normalized_DESEQ2.logtransformed.txt")
genes <- scan("Genes-for-Enrichment_T.leiboldiana_unsigned8_vst10cs4_yellowgreen.txt", character(), quote = "")
geneInfo <- read.table("geneInfo_co-expression_Tfasc_vsd_10c4s_signed18.csv")
goterms <- read.delim("GO-term_enrichment_T.leiboldiana_unsigned8_mod-yellowgreen.txt", sep = "\t", header = T)

module_counts<- subset(counts, rownames(counts) %in% genes)
setdiff(genes, rownames(module_counts))

mean_count <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count) <- c("0100", "0500", "0900","1300", "1700", "2100")
mean_count$gene_id <- row.names(mean_count)
mean_count_m <- melt(mean_count, id.vars = "gene_id")
colnames(mean_count_m) <- c("gene_id", "time", "count")

module <- str_split("Genes-for-Enrichment_T.fasciculata_signed18_vst10c4s_lightcyan.txt", "\\_|\\.")[[1]][6]
species <- str_split("Genes-for-Enrichment_T.fasciculata_signed18_vst10c4s_lightcyan.txt", "\\_|\\.")[[1]][3]
network <- str_split("Genes-for-Enrichment_T.fasciculata_signed18_vst10c4s_lightcyan.txt", "\\_|\\.")[[1]][4] 

module_geneInfo <- subset(geneInfo, rownames(geneInfo) %in% genes)
membership <- as.data.frame(cbind((rownames(module_geneInfo)), 
                                  (module_geneInfo$MM.lightcyan)))
membership$V2 <- as.numeric(membership$V2)
summary(membership)
colnames(membership) <- c("gene_id", "membership")

malate <- goterms[grep("malate", goterms$Term), 6]
malate <- unlist(str_split(malate, ", "))
malate_expr <- subset(mean_count_m, gene_id %in% malate)

phospho <- goterms[grep("phosphoenolpyruvate", goterms$Term), 6]
phospho <- unlist(str_split(phospho, ", "))
phospho_expr <- subset(mean_count_m, gene_id %in% phospho)

vacuole <- goterms[grep("vacuolar acidification", goterms$Term), 6]
vacuole <- unlist(str_split(vacuole, ", "))
vacuole_expr <- subset(mean_count_m, gene_id %in% vacuole)

stomata <- goterms[grep("stomata", goterms$Term), 6]
stomata <- unlist(str_split(stomata, ", "))
stomata_expr <- subset(mean_count_m, gene_id %in% stomata)

mm <- merge(mean_count_m, membership, by = "gene_id")
mean_count_m <- mm
mean_count_m$membership <- as.numeric(mean_count_m$membership)

means <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"), 
               count = c(mean(mean_count$`0100`), mean(mean_count$`0500`),
                 mean(mean_count$`0900`), mean(mean_count$`1300`),
                 mean(mean_count$`1700`), mean(mean_count$`2100`)))


ggplot(mean_count_m, aes(x=time, y=count)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.3)

# Colour in grey
ggplot(mean_count_m, aes(x=time, y=count, group = gene_id)) + 
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means, aes(group = 1), size = 1.25, color = "black")

# Colour by degree of membership
ggplot(mean_count_m, aes(x=time, y=count, group = gene_id, color = membership)) + 
  geom_point() +
  geom_line() +
  geom_line(data = means, aes(group = 1), size = 1.25, color = "black") +
  scale_color_gradient(low="blue", high="red") +
  ylim(0,15000)

# Highlight genes of interest
ggplot(mean_count_m, aes(x=time, y=count, group = gene_id)) + 
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means, aes(group = 1), size = 1, color = "black") +
  geom_line(data = malate_expr, aes(group = gene_id), size = 1, color = "darkgreen")

pdf(paste("Expression_curve_", species, "_", network, "_", module,".pdf"), width = 8, height = 10)
cols <- c('Malate'='red', 'PhosphoEnolPyruvate'='blue', 'Vacuole'='purple', 'Stomata'='darkgreen')
ggplot(mean_count_m, aes(x=time, y=count, group = gene_id)) + 
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means, aes(group = 1), size = 1, color = "black") +
  geom_line(data = malate_expr, aes(group = gene_id, color = "Malate"), size = 1) +
  geom_line(data = phospho_expr, aes(group = gene_id), size = 1, color = "blue") +
  geom_line(data = vacuole_expr, aes(group = gene_id), size = 1, color = "purple") + 
  geom_line(data = stomata_expr, aes(group = gene_id), size = 1, color = "darkgreen") +
  ggtitle(paste("Expression curve of ",module, "in", species, "_", network)) +
  scale_colour_manual(name='Function of Interest',
                     values=cols)+
  theme(legend.title=element_text(size=10
                                  ),
        legend.text=element_text(size=10))
  dev.off()  

       