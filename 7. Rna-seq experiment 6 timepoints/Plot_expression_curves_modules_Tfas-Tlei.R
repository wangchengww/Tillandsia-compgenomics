setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/')
library("ggplot2")
library("reshape2")
library("stringr")
library(grid)
library(gridExtra)

counts <- read.table("counts.Tfas_Tlei_6_timepoints.filtr-normalized_DESEQ2.txt", header = T, row.names = 1)
genes <- scan("Genes-for-Enrichment_Tfas-Tlei_unsigned7_vst_violet.txt", character(), quote = "")
geneInfo <- read.table("geneInfo_co-expression_Tfasc_vsd_10c4s_signed18.csv")
goterms <- read.table("GO-term_enrichment_Tfas-Tlei_unsigned7_mod-violet.txt", sep = "\t", header = T)

module_counts<- subset(counts, rownames(counts) %in% genes)
module_counts_Tfas <- module_counts[, c(1:36)]
module_counts_Tfas <- module_counts_Tfas[, -c(12,36)]
module_counts_Tlei <- module_counts[, c(37:72)]

mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tfas) <- c("0100", "0500", "0900","1300", "1700", "2100")
mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tlei) <- c("0100", "0500", "0900","1300", "1700", "2100")

mean_count_Tfas$gene_id <- row.names(mean_count_Tfas)
mean_count_Tlei$gene_id <- row.names(mean_count_Tlei)

summary(mean_count_Tfas)

mean_count_Tfas_m <- melt(mean_count_Tfas, id.vars = "gene_id")
colnames(mean_count_Tfas_m) <- c("gene_id", "time", "count")
mean_count_Tlei_m <- melt(mean_count_Tlei, id.vars = "gene_id")
colnames(mean_count_Tlei_m) <- c("gene_id", "time", "count")

module_geneInfo <- subset(geneInfo, rownames(geneInfo) %in% genes)
membership <- as.data.frame(cbind((rownames(module_geneInfo)),
                                  (module_geneInfo$MM.lightcyan)))
membership$V2 <- as.numeric(membership$V2)
summary(membership)
colnames(membership) <- c("gene_id", "membership")

malate <- goterms[grep("malate", goterms$Term), 6]
malate <- unlist(str_split(malate, ", "))
malate_expr_Tfas <- subset(mean_count_Tfas_m, gene_id %in% malate)
malate_expr_Tlei <- subset(mean_count_Tlei_m, gene_id %in% malate)


mm <- merge(mean_count_m, membership, by = "gene_id")
mean_count_m <- mm
mean_count_m$membership <- as.numeric(mean_count_m$membership)

means_Tfas <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
               count = c(mean(mean_count_Tfas$`0100`), mean(mean_count_Tfas$`0500`),
                 mean(mean_count_Tfas$`0900`), mean(mean_count_Tfas$`1300`),
                 mean(mean_count_Tfas$`1700`), mean(mean_count_Tfas$`2100`)))

means_Tlei <- data.frame(time=c("0100", "0500","0900", "1300", "1700", "2100"),
                         count = c(mean(mean_count_Tlei$`0100`), mean(mean_count_Tlei$`0500`),
                                   mean(mean_count_Tlei$`0900`), mean(mean_count_Tlei$`1300`),
                                   mean(mean_count_Tlei$`1700`), mean(mean_count_Tlei$`2100`)))
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
p1 <- ggplot(mean_count_Tfas_m, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means_Tfas, aes(group = 1), size = 1, color = "black") +
  geom_line(data = malate_expr_Tfas, aes(group = gene_id), size = 1, color = "darkgreen") +
  scale_y_continuous(name="",
                     breaks = c(seq(from = 0, to = 25000, by = 5000)),
                     limits = c(0,25000))

p2 <- ggplot(mean_count_Tlei_m, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means_Tlei, aes(group = 1), size = 1, color = "black") +
  geom_line(data = malate_expr_Tlei, aes(group = gene_id), size = 1, color = "darkgreen") +
  scale_y_continuous(name="",
                     breaks = c(seq(from = 0, to = 25000, by = 5000)),
                     limits = c(0,25000))

grid.arrange(p1, p2, nrow = 1, top = "Expression curve of MEViolet in T. fasciculata
             (left) and T. leiboldiana (right)")
