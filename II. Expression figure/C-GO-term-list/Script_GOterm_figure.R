setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/C-GO-term-list/')
setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/II. Expression figure/C-GO-term-list/')

library("GOplot")
library(stringr)
go <- read.delim("GO-term_enrichment_mod_TLEI-REF_exonic.txt", header = T, sep = "\t")
ortho_info <- read.delim("orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.blastandsearch_noAcom.txt",
                          sep = "\t")
colnames(ortho_info) <- c("gene_ID", "chr", "start", "end", "GOterm", "Description", "orthogroup", "count_Acom", 
                          "count_Tfas", "count_Tlei")
genelist <- unlist(str_split(go$Genes, ", "))
genes <- subset(ortho_info, gene_ID %in% genelist)

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

GOBar <- function(data, display, order.by.zscore = T, title, zsc.col){
  id <- adj_pval <- zscore <- NULL
  if (missing(display)) display <- 'single'
  if (missing(title)) title <- ''
  if (missing(zsc.col)) zsc.col <- c('firebrick1', 'white', 'dodgerblue1')
  colnames(data) <- tolower(colnames(data))
  data$adj_pval <- -log(data$adj_pval, 10)
  sub <- data[!duplicated(data$term), ]
  if (order.by.zscore == T) {
    sub <- sub[order(sub$zscore, decreasing = T), ]
    g <-  ggplot(sub, aes(y = factor(term, levels = rev(stats::reorder(term, adj_pval))), x = adj_pval, fill = zscore)) +
      geom_bar(stat = 'identity', colour = 'black') +
      scale_fill_gradient2('Per-family\ngene count\n(z-score)\n', space = 'Lab', low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], guide = guide_colourbar(title.position = "top", title.hjust = 0), 
                           breaks = c(min(sub$zscore), max(sub$zscore)), labels = c('T.fas = T. lei', 'T.fas > T. lei')) +
      labs(title = title, y = '', x = '-log (adj p-value)') +
      geom_text(aes(y=term, x=-0.25, label=count))
  }else{
    sub <- sub[order(sub$adj_pval, decreasing = T), ]
    leg <- theme(legend.justification = c(1, 1), legend.position = c(2, 2.1), legend.background = element_rect(fill = 'transparent'),
                 legend.box = 'vertical', legend.direction = 'horizontal')
    g <-  ggplot(sub, aes( x = factor(id, levels = reorder(id, adj_pval)), y = zscore, fill = adj_pval)) +
      geom_bar(stat = 'identity', colour = 'black') +
      scale_fill_gradient2('Significance', space = 'Lab', low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], guide = guide_colourbar(title.position = "top", title.hjust = 0.5), breaks = c(min(sub$adj_pval), max(sub$adj_pval)), labels = c('low', 'high')) +
      labs(title = title, x = '', y = 'z-score') +
      leg
  }
  if (display == 'single'){
    g + theme(axis.text.x = element_text(angle = 360, vjust = 0.5, hjust=1), axis.line = element_line(colour = 'grey80'), axis.ticks = element_line(colour = 'grey80'),
              axis.title = element_text(size = 12, face = 'bold'), axis.text = element_text(size = 12), panel.background = element_blank(), 
              panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())        
  }else{
    g + facet_grid(.~category, space = 'free_x', scales = 'free_x')+
      theme(axis.text.x = element_text(angle = 180, vjust = 0.5, hjust=1), axis.line = element_line(colour = 'grey80'), axis.ticks = element_line(colour = 'grey80'),
            axis.title = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 14), panel.background = element_blank(), 
            panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  }
}
genes <- genes[,c(1,9,10)]
genes$difference <- genes$count_Tfas - genes$count_Tlei
genes <- genes[,c(1,4)]
colnames(genes) <- c("ID", "logFC")
circ <- circle_dat(go, genes)
circ_sub <- circ[grepl("oxaloacetate", circ$term) | grepl("glycoly", circ$term) | grepl("circadian", circ$term)
                 | grepl("fructo", circ$term) | grepl("malate", circ$term) | grepl("starch", circ$term) 
                 | grepl("vacuol", circ$term)| grepl("tricarboxylic", circ$term) | grepl("heat", circ$term)
                 | grepl("ATPase", circ$term),]
pdf("GOterm_CAM_list.pdf", width = 12, height = 9)
GOBar(circ_sub, zsc.col = c("#1e6091", "#52b69a", "#d9ed92"))
dev.off()
GOBubble(circ, labels = 3)


ggplot(circ_sub, aes(y=term, x=-log(adj_pval,10))) +
  # add bars
  geom_col(aes(fill=zscore), colour="black") +
  
  # colour bars with gradient
  scale_fill_gradient2('z-score \nper-family gene count\n',low="#d9ed92", mid = "#52b69a", 
                       high="#1e6091",labels=c('T.fas = T. lei', 'T.fas > T. lei'),
                       breaks=c(min(circ_sub$zsc), max(circ_sub$zsc)), limits=c(min(circ_sub$zsc), max(circ_sub$zsc))) +
  # add Count of genes
  geom_text(aes(y=term, x=-0.25, label=count))
  
  # add grouping of GO terms
  geom_point(data=d[!(is.na(d$GO_cl)),], aes(y=Term, x=-0.6, colour=GO_cl, shape=GO_cl), size=4) +
  
  # change shape and colour of GO term groupings and merge into one legend
  scale_shape_manual(name="Grouping of GO terms", values=rep(c(15,16,17), each=3)) +
  scale_colour_manual(name = "Grouping of GO terms", values = rep(c("grey65", "gold", "black"),3))

# modify area to place labels outside axis lines 
coord_cartesian(clip = 'off', ylim = c(0, nrow(d)+1), xlim=c(0, 5.7)) +
  
  # change space around the plot area
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand=c(t=0, r=0.75, b=0, l=0)) +
  
  # set the theme
  theme(plot.margin=unit(c(t=10, r=10, b=10, l=1),"mm"),
        panel.background=element_rect(colour="white", fill="white", size=0.5),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.ticks.y=element_blank(), axis.title=element_blank(),
        axis.text.y=element_text(size=14), axis.text.x=element_text(size=14),
        legend.text=element_text(size=14), legend.title=element_text(size=16, face="bold"),
        legend.key=element_blank(), legend.margin=margin(unit(c(t=40, r=1, b=1, l=1),"mm")),
        legend.spacing = unit(5, "mm"))