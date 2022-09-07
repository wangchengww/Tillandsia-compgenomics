setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/C-GO-term-list/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/C-GO-term-list/')

library("GOplot")
library(stringr)

### Define functions
###########################################################################

# Create input dataset
circle_dat <- function(terms, genes){
  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), ', ')
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), ',')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
  if(class(logFC) == 'factor'){
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1; zsc <- c()
  for (c in 1:length(count)){
    value <- 0
    e <- s + count[c] - 1
    value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, ifelse(x == 0, 0, -1)))
    zsc <- c(zsc, sum(value) / sqrt(count[c]))
    s <- e + 1
  }
  if (is.null(terms$id)){
    df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}

# P-value is on x-axis and z-score is colour scale
GOBar <- function(data, display, order.by.zscore = T, title, zsc.col){
  id <- adj_pval <- zscore <- NULL
  if (missing(zsc.col)) zsc.col <- c('firebrick1', 'white', 'dodgerblue1')
  colnames(data) <- tolower(colnames(data))
  data$adj_pval <- -log(data$adj_pval, 10)
  sub <- data[!duplicated(data$term), ]
  if (order.by.zscore == T) {
    sub <- sub[order(sub$zscore, decreasing = T), ]
    g <-  ggplot(sub, aes(y = factor(term, levels = rev(stats::reorder(term, adj_pval))), x = adj_pval, fill = zscore)) +
      geom_bar(stat = 'identity', colour = 'black') +
      theme_bw()+
      scale_fill_gradient2('Per-family\ngene count\n(z-score)\n', space = 'Lab', low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], guide = guide_colourbar(title.position = "top", title.hjust = 0), 
                           breaks = c(min(sub$zscore), max(sub$zscore)), labels = c('T.fas < T. lei', 'T.fas > T. lei')) +
      labs(title = '', y = '', x = '-log (adj p-value)') +
      geom_text(aes(y=term, x=-0.25, label=count),size = 4) +
      theme(axis.text=element_text(size=10))
    g
  }
}

# Z-score is on x-axis and p-value is colour scale
GOBar <- function(data, display, order.by.zscore = T, title, zsc.col){
  id <- adj_pval <- zscore <- NULL
  if (missing(zsc.col)) zsc.col <- c('firebrick1', 'white', 'dodgerblue1')
  colnames(data) <- tolower(colnames(data))
  data$adj_pval <- -log(data$adj_pval, 10)
  sub <- data[!duplicated(data$term), ]
  if (order.by.zscore == T) {
    sub <- sub[order(sub$zscore, decreasing = T), ]
    g <-  ggplot(sub, aes(y = factor(term, levels = rev(stats::reorder(term, adj_pval))), x = zscore, fill = adj_pval)) +
      geom_bar(stat = 'identity', colour = 'black') +
      theme_bw()+
      scale_fill_gradient2('Per-family\ngene count\n(z-score)\n', space = 'Lab', low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], guide = guide_colourbar(title.position = "top", title.hjust = 0), 
                           breaks = c(min(sub$adj_pval), max(sub$adj_pval)), labels = c('T.fas < T. lei', 'T.fas > T. lei')) +
      labs(title = '', y = '', x = 'z-score') +
      geom_text(aes(y=term, x=-0.25, label=count),size = 4) +
      theme(axis.text=element_text(size=10))
    g
  }
}

##########################################################################

### Read data
go <- read.delim("GO-term_enrichment.COMBINED.Tfas-Tlei.onepvalue.txt", header = T, sep = "\t")
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

### Create input data
circ <- circle_dat(go, genes)
circ_sub <- circ[(grepl("oxaloacetate", circ$term) | grepl("glycoly", circ$term) | grepl("circadian", circ$term)
                  | grepl("fructo", circ$term) | grepl("malate", circ$term) | grepl("starch", circ$term) 
                  | grepl("vacuol", circ$term)| grepl("tricarboxylic", circ$term) | grepl("heat", circ$term)
                  | grepl("ATPase", circ$term)| grepl("salt", circ$term)| grepl("photoperiodism", circ$term)| 
                    grepl("osmotic", circ$term)) & (circ$category == "BP" | circ$category == "MF"),]

### Plot figure
pdf("GOterm_CAM_list_Tfas-Tlei_COMBINED.onepvalue2.pdf", width = 8, height = 12)
GOBar(circ_sub, zsc.col = c("#1e6091", "#52b69a", "#d9ed92"))
dev.off()
