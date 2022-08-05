### Assessment of gene duplications in Tfas ###

### Here I check whether inferred gene duplication in Tillandsia fasciculata are legitimate or not
### based on median coverage calculations after mapping raw illumina data, 50x, to the reference of
### Tfas. The idea is that genuine duplications will maintain stable coverage which is similar to
### the background, whereas false duplications will have lower coverage.

library(ggplot2)
library(plyr)
library(dplyr)
library(wesanderson)
library(reshape2)

setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/5")
setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/5. Gene-family-evo-tfas-tlei/")

tfas_genes <- read.table("Tfas_pergene_mediancov_and_orthoinfo.txt", header = T)

mean(tfas_genes$median_cov) # 47.5
mean(tfas_genes$mean_cov) # 48.04

# Set colours categories
# mycolors <- c(wes_palette("Darjeeling1")[1], wes_palette("Darjeeling1")[2], wes_palette("FantasticFox1")[3],
#               wes_palette("Darjeeling1")[3], wes_palette("Rushmore1")[3], wes_palette("GrandBudapest1")[3])
# Divide genes into categories based on copy number
#dup <- c()
#for (i in 1:nrow(tfas_genes)){
#  if (tfas_genes[i,6] == 1){
#    dup <- c(dup, "single-copy")
#  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,5] == 0 && tfas_genes[i,7] == 0){
#    dup <- c(dup, "unique-multi-copy")
#  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,6] > tfas_genes[i,5] && tfas_genes[i,6] > tfas_genes[i,7]){
#    dup <- c(dup, "largest-multi-copy")
#  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,6] < tfas_genes[i,5] && tfas_genes[i,6] < tfas_genes[i,7]){
#    dup <- c(dup, "smallest-multi-copy")
#  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,6] == tfas_genes[i,5] && tfas_genes[i,6] == tfas_genes[i,7]){
#    dup <- c(dup, "ancestral-multi-copy")
#  } else {
#    dup <- c(dup, "middle-multi-copy")
#  }
#}
#tfas_genes$duplicated <- dup

# simplified categories
dup_tfas <- c()
for (i in 1:nrow(tfas_genes)){
  if (tfas_genes[i,6] == 1 && tfas_genes[i,5] == 1 && tfas_genes[i,7] ==1){
    dup_tfas <- c(dup_tfas, "ancestral-single-copy")
  } else if (tfas_genes[i,6] == 1){
    dup_tfas <- c(dup_tfas, "single-copy")
  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,5] == 0 && tfas_genes[i,7] == 0){
    dup_tfas <- c(dup_tfas, "unique-multi-copy")
  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,6] == tfas_genes[i,5] && tfas_genes[i,6] == tfas_genes[i,7]){
    dup_tfas <- c(dup_tfas, "ancestral-multi-copy")
  } else {
    dup_tfas <- c(dup_tfas, "multi-copy")
  }
}
tfas_genes$duplicated <- dup_tfas

# Obtain average coverage for each category
mu <- ddply(tfas_genes, "dup", summarise, grp.mean=mean(median_cov))
head(mu)
# dup grp.mean
# ancestral-multi-copy    40.28850
# ancestral-single-copy   42.91895
# multi-copy              42.65542
# single-copy             65.12410
# unique-multi-copy       52.45243

#Statistics for the full genome
# count    6.002822e+08
# mean     46.17124
# std      110.3986
# min      0
# 25%      21.00000
# 50%      36.00000
# 75%      51.00000
# max      8070.000

# Set colors for simple categories
mycolors <- c(wes_palette("Darjeeling1")[1], wes_palette("Darjeeling1")[2], wes_palette("FantasticFox1")[3],
              wes_palette("Darjeeling1")[3], wes_palette("GrandBudapest1")[3])

# Density plot of mean coverage per gene per category
ggplot(tfas_genes, aes(x=mean_cov, color=dup)) +
  geom_density() +
  xlim(0,150) +
  geom_vline(xintercept = 46.1712, linetype = "longdash", colour = "gray28") +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)

# Density plot of median coverage per gene per category
ggplot(tfas_genes, aes(x=median_cov, color=dup)) +
  geom_density() +
  xlim(0,150) +
  geom_vline(xintercept = 46.1712, linetype = "longdash", colour = "gray28") +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)

---------------
# UPDATE 27.07: I skipped establishing a cutoff and corrected sizes for all multicopy families
#
### Establish a cutoff to separate "real" from "faulty" gene models
## Using the package "cutoff", we can establish a cut-off value that separates the bimodal
## distribution into two unimodal, normal distributions. The cut-off value represents the point
## where the likelihood of a datapoint belonging only to one of the two distributions is high
## enough that we can call it with confidence. This will be a slightly conservative value, but
## allows us to separate the two groups based on statistics and probability.

#devtools::install_github("choisy/cutoff")
#install.packages("cli")
#install.packages("bbmle")
#library(cutoff)

#tfas_multicopy_mean_cov <- as.vector(tfas_genes[tfas_genes$duplicated == 'multi-copy', 3])
#tfas_multicopy_mean_cov_df <- tfas_genes[tfas_genes$duplicated == 'multi-copy',]
#tfas_multicopy_mean_cov_subset <- tfas_multicopy_mean_cov[tfas_multicopy_mean_cov < 100]
#tfas_multicopy_FMM <- em(tfas_multicopy_mean_cov_subset, 'normal', 'normal')
#hist(tfas_multicopy_mean_cov_subset,10000,F,xlim=c(0,100),xlab="median coverage",ylab="density",main=NULL,col="grey")
#lines(tfas_multicopy_FMM,lwd=1.5,col="red")
#abline(v=34.5,lty=2,col="lightblue")

#cut_off <- cutoff(tfas_multicopy_FMM)
#cut_off # Genes with mean cov < 34.5 belong to "faulty" category

#faulty_genes <- tfas_multicopy_mean_cov_df[tfas_multicopy_mean_cov_df$mean_cov < 34.5,] # 4753 genes
#table(faulty_genes$Tfas_count)

# Recording all genes with high coverage to explore further (mostly mitochondrial and plastid)
# high_cov_genes <- tfas_genes[tfas_genes$mean_cov > 100,]
# write.table(high_cov_genes, file = "100x_genes_Tfas.txt")

# Recording all genes with a mean coverage < 34.5 to correct gene family sizes
#write.table(faulty_genes, file = "faulty_gene_families_35x.txt")

# Select all orthogroups containing faulty genes
# library(dplyr)
# orthogroups_w_faulty_genes <- tfas_genes %>%
#  filter(og_id %in% faulty_genes$og_id)
# This resulted in 6861 genes and 1981 orthogroups
-----------------
  
# Correcting the gene family sizes based on coverage differences
  
# Now I select just multicopy genes - this immediately captures the full orthogroup as well.
# At this step, I remove all genes that come from orthogroups containing plastid, mitochondrial or ribosomal annotations.
# The reasoning behind this is that we can't truly use coverage statistics to adjust sizes since plastid and mitochondrial
# DNA is usually sequenced at high coverage (there are thousands of copies in each cell). Additionally, we will skip 
# looking into gene family evolution of these types of genes since they have quite different dynamics (especially 
# ribosomal genes).
tfas_multicopy_genes <- tfas_genes[tfas_genes$duplicated == 'multi-copy',] # 8400 genes, 2632 orthogroups
mito_og <- read.table("mito_plastid_ribo_OGs_to_remove.txt", header = F, sep = "\t") # list with mito-ribo-plastid orthogroups
colnames(mito_og) <- c("og_id")
tfas_multicopy_genes_no_mito <- tfas_multicopy_genes %>%
  filter(!(og_id %in% mito_og$og_id)) # 8139 genes, 2556 orthogroups

# Correction without accounting for coverage variability
corr_Tfas_simple <- data.frame()
for (i in unique(tfas_multicopy_genes_no_mito$og_id)){
  orthogroup <- tfas_multicopy_genes_no_mito[tfas_multicopy_genes_no_mito$og_id == i,]
  nr_genes <- as.integer(orthogroup[1,6])
  total_mean_cov <- sum(orthogroup$mean_cov)
  expected_mean_cov <- nr_genes * 42.9189
  correction_factor <- total_mean_cov/expected_mean_cov
  new_size = round((nr_genes * correction_factor), digits = 0)
  #if (correction_factor > 1){
  #  new_size = nr_genes
  #}
  if (total_mean_cov < 42.9189){
    new_size = 1
  }
  corr_family_size <- cbind(i, as.numeric(correction_factor), as.integer(new_size),
                            as.integer(nr_genes), as.numeric(total_mean_cov), expected_mean_cov)
  corr_Tfas_simple <- rbind(corr_Tfas_simple, corr_family_size)
}
colnames(corr_Tfas_simple) <- c("og_id", "correction_factor", "corr_Tfas_count","old_Tfas_count",
                         "total_mean_cov", "expected_mean_cov")

# Correction accounting for variability
quantile(tfas_genes[tfas_genes$duplicated == "ancestral-single-copy", 3], c(.01, .05, .1, .25, .5,.75, .9, .95, .99)) 
#1%       5%      10%      25%      50%      75%      90%      95%      99% 
#15.13276 22.74074 27.60464 34.74492 39.95281 43.76340 46.56272 48.43682 69.65468 
lower_thresh=34.74492# 25th quantile of the average coverage of ancestral single copy genes to account for variation
upper_thresh=43.76340 # 75th quantile of the average coverage of ancestral single copy genes
corr_Tfas_var <- data.frame()
for (i in unique(tfas_multicopy_genes_no_mito$og_id)){
  orthogroup <- tfas_multicopy_genes_no_mito[tfas_multicopy_genes_no_mito$og_id == i,]
  nr_genes <- as.integer(orthogroup[1,6])
  total_mean_cov <- sum(orthogroup$mean_cov)
  lower_thresh_total <- nr_genes * lower_thresh
  upper_thresh_total <- nr_genes * upper_thresh
  expected_mean_cov <- nr_genes * 42.9189
  correction_factor <- total_mean_cov/expected_mean_cov
  new_size = round((nr_genes * correction_factor), digits = 0)
  if (lower_thresh_total < total_mean_cov && total_mean_cov < upper_thresh_total){
    new_size = nr_genes
  }
  #if (correction_factor > 1){
  #  new_size = nr_genes
  #}
  if (total_mean_cov < 42.9189){
    new_size = 1
  }
  corr_family_size <- cbind(i, as.numeric(correction_factor), as.integer(new_size),
                            as.integer(nr_genes), as.numeric(total_mean_cov), expected_mean_cov, 
                            lower_thresh_total,upper_thresh_total)
  corr_Tfas_var <- rbind(corr_Tfas_var, corr_family_size)
}
colnames(corr_Tfas_var) <- c("og_id", "correction_factor", "corr_Tfas_count_PO","old_Tfas_count_PO", "total_mean_cov", "expected_mean_cov", "Lower_thresh", "Upper_thresh")

nrow(corr_Tfas_var[(corr_Tfas_var$total_mean_cov > corr_Tfas_var$Lower_thresh) & (corr_Tfas_var$total_mean_cov < corr_Tfas_var$Upper_thresh),])
# 549 orthogroups are inside threshold 
nrow(corr_Tfas_var[corr_Tfas_var$corr_Tfas_count_PO == corr_Tfas_var$old_Tfas_count_PO,]) 
# 825 genes have been unchanged (so more than genes in thresholds, probably because of rounding) 
nrow(corr_Tfas_simple[corr_Tfas_simple$corr_Tfas_count == corr_Tfas_simple$old_Tfas_count,])
# 764 orthogroups were unchanged when not accounting for variability, so accounting for variability only affects 61 
# orthogroups or 2.6 % of the dataset.

# Per-gene size correction: I tested a system where size corrections is done on a per-gene basis. 
# I assign a probablity of each gene based on its coverage with respect to the mean coverage of 
# ancestral single copy genes. If the gene has the same coverage, it will have a prob of 1 and counted as
# a full gene. Lower coverages lead to adusted weighting. At the end, the adjusted weights are summed up
# into a new orthogroup size.
corr_Tfas_per_gene <- data.frame()
for (i in unique(tfas_multicopy_genes_no_mito$og_id)){
  orthogroup <- tfas_multicopy_genes_no_mito[tfas_multicopy_genes_no_mito$og_id == i,]
  weighting_vect <- c()
  for (j in 1:nrow(orthogroup)){
    mean_cov <- orthogroup[j,3]
    weighting <- mean_cov / 42.9189
    weighting_vect <- c(weighting_vect, weighting)
  }
  new_size = round(sum(weighting_vect))
  nr_genes <- as.integer(orthogroup[1,6])
  corr_family_size <- cbind(i, as.integer(new_size),
                            as.integer(nr_genes))
  corr_Tfas_per_gene <- rbind(corr_Tfas_per_gene, corr_family_size)
}
colnames(corr_Tfas_per_gene) <- c("og_id", "corr_Tfas_count_PG","old_Tfas_count_PG")

# Comparing a per-gene correction approach to a per-orthogroup approach. Here, corrections were made
# in both directions (increase and decrease in size) and without accounting for variability
comparison_siye_corr <- merge(corr_Tfas_simple, corr_Tfas_per_gene, by = "og_id")
comparison_size_corr <- comparison_siye_corr[,c(1,4,3,7)]
colnames(comparison_size_corr) <- c("og_id", "uncorrected_count", "correction_per_OG", "correction_per_GENE")
# calculate differences in corrections
comparison_size_corr$diff_per_OG <- as.numeric(comparison_size_corr$uncorrected_count) - as.numeric(comparison_size_corr$correction_per_OG) 
comparison_size_corr$diff_per_GENE <- as.numeric(comparison_size_corr$uncorrected_count) - as.numeric(comparison_size_corr$correction_per_GENE) 
comparison_size_corr$diff_approach <- as.numeric(abs(comparison_size_corr$diff_per_OG)) - abs(as.numeric(comparison_size_corr$diff_per_GENE))
# See how much the differences deviate 
diff_melt <- melt(comparison_size_corr[, c(1,5,6)], id.vars = "og_id")
ggplot(diff_melt, aes(x=value, color=variable)) +
  geom_histogram(fill = "white", alpha=0.5, position="identity", binwidth = 1) +
  xlim(-25,25) + 
  xlab("Differences in size corrections (difference_per_orthogroup - difference_per_gene")

ggplot(comparison_size_corr, aes(x=diff_approach)) +
  geom_boxplot()
mean(comparison_size_corr$diff_approach) # -0.013
table(comparison_size_corr$diff_approach)

# The per-gene and per-orthogroup approaches lead to very similar size corrections. 2323 orthogroups were corrected in 
# the same way, 30 had a -1 difference between the two approaches.

# Given the little differences between per-orthogroup and per-gene method, and the little importance 
# and limitations of accounting for variability in coverage, I decided to stick to the corrected sizes of per-orthogroup
# approach without accounting for variability.
write.table(corr_Tfas_simple, file = "corrected_family_sizes_Tfas.txt", sep = "\t", 
            quote = F, row.names = F)

###########################################################################################
### Same for Tleiboldiana
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
tlei_genes <- read.table("Tlei_pergene_mediancov_and_orthoinfo.txt", header = T)
mean(tlei_genes$median_cov) # 67.04
mean(tlei_genes$mean_cov) # 66.73

# dup <- c()
# for (i in 1:nrow(tlei_genes)){
#   if (tlei_genes[i,7] == 1){
#     dup <- c(dup, "single-copy")
#   } else if (tlei_genes[i,7] > 1 && tlei_genes[i,5] == 0 && tlei_genes[i,6] == 0){
#     dup <- c(dup, "unique-multi-copy")
#   } else if (tlei_genes[i,7] > 1 && tlei_genes[i,7] > tlei_genes[i,5] && tlei_genes[i,7] > tlei_genes[i,6]){
#     dup <- c(dup, "largest-multi-copy")
#   } else if (tlei_genes[i,7] > 1 && tlei_genes[i,7] < tlei_genes[i,5] && tlei_genes[i,7] < tlei_genes[i,6]){
#     dup <- c(dup, "smallest-multi-copy")
#   } else if (tlei_genes[i,7] > 1 && tlei_genes[i,7] == tlei_genes[i,5] && tlei_genes[i,7] == tlei_genes[i,6]){
#     dup <- c(dup, "ancestral-multi-copy")
#   } else {
#     dup <- c(dup, "middle-multi-copy")
#   }
# }
# tlei_genes$duplicated <- dup

# simplified
dup_tlei <- c()
for (i in 1:nrow(tlei_genes)){
  if (tlei_genes[i,7] == 1 && tlei_genes[i,5] == 1 && tlei_genes[i,6] ==1){
    dup_tlei <- c(dup_tlei, "ancestral-single-copy")
  } else if (tlei_genes[i,7] == 1){
    dup_tlei <- c(dup_tlei, "single-copy")
  } else if (tlei_genes[i,7] > 1 && tlei_genes[i,5] == 0 && tlei_genes[i,6] == 0){
    dup_tlei <- c(dup_tlei, "unique-multi-copy")
  } else if (tlei_genes[i,7] > 1 && tlei_genes[i,7] == tlei_genes[i,5] && tlei_genes[i,7] == tlei_genes[i,6]){
    dup_tlei <- c(dup_tlei, "ancestral-multi-copy")
  } else {
    dup_tlei <- c(dup_tlei, "multi-copy")
  }
}
tlei_genes$duplicated <- dup_tlei

mu <- ddply(tlei_genes, "dup", summarise, grp.mean=mean(mean_cov))
head(mu)
# dup grp.mean
# 1  ancestral-multi-copy 49.28971
# 2 ancestral-single-copy 50.69818
# 3            multi-copy 79.84922
# 4           single-copy 91.58547
# 5     unique-multi-copy 58.87240

ggplot(tlei_genes, aes(x=mean_cov, colour=duplicated)) +
  geom_density() +
  xlim(0,150) +
  geom_vline(xintercept = 53.34, linetype = "longdash", color = "gray28") +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)

ggplot(tlei_genes, aes(x=median_cov, color=dup_tlei)) +
  geom_density() +
  xlim(0,150) +
  geom_vline(xintercept = 53.34, linetype = "longdash", colour = "gray28") +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)

# Corrections on all multicopy genes, except for ribosomal / plastid / mitochondrial.
tlei_multicopy_genes <- tlei_genes[tlei_genes$duplicated == 'multi-copy',] # 4306 genes, 1450 orthogroups
tlei_multicopy_genes_no_mito <- tlei_multicopy_genes %>%
  filter(!(og_id %in% mito_og$og_id)) # 4098 genes, 1377 orthogroups

corr_Tlei <- data.frame()
for (i in unique(tlei_multicopy_genes_no_mito$og_id)){
  orthogroup <- tlei_multicopy_genes_no_mito[tlei_multicopy_genes_no_mito$og_id == i,]
  nr_genes <- as.integer(orthogroup[1,7])
  total_mean_cov <- sum(orthogroup$mean_cov)
  expected_mean_cov <- nr_genes * 50.69818
  correction_factor <- total_mean_cov/expected_mean_cov
  new_size = round((nr_genes * correction_factor), digits = 0)
  # if (correction_factor > 1){
  #   new_size = nr_genes
  # }
  if (total_mean_cov < 50.69818){
    new_size = 1
  }
  corr_family_size <- cbind(i, as.numeric(correction_factor), as.integer(new_size),
                            as.integer(nr_genes), as.numeric(total_mean_cov), expected_mean_cov)
  corr_Tlei <- rbind(corr_Tlei, corr_family_size)
}
colnames(corr_Tlei) <- c("og_id", "correction_factor", "corr_Tlei_count","old_Tlei_count", "total_mean_cov", "expected_mean_cov")

write.table(corr_Tlei, file = "corrected_family_sizes_Tlei.txt", sep = "\t", 
            quote = F, row.names = F)

### Making supplementary figure ###

pacman::p_load("ggplot2","grid", "gridExtra","RColorBrewer")

# Extract legend to plot legend for multiple plots
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract x-axis
titleGrob <- function(label, margin=unit(c(b=1, l=1, t=0.1, r=1), "line"), ..., debug=FALSE){
  library(gtable)
  tg <- textGrob(label, ...)
  w <- grobWidth(tg)
  h <- grobHeight(tg)
  
  g <- gtable("title",
              widths = unit.c(margin[2], w, margin[4]), 
              heights = unit.c(margin[3], h, margin[1]), respect = FALSE)
  if(debug){
    rg <- rectGrob()
    pos <- expand.grid(x=1:3, y=1:3)[-5,]
    g <- gtable_add_grob(g, list(rg, rg, rg, rg, rg, rg, rg, rg), t=pos$y, l=pos$x)
  }
  gtable_add_grob(g, tg, t=2, l = 2)
}

p1 <- ggplot(tfas_genes, aes(x=mean_cov, color=dup_tfas)) +
  geom_density() +
  theme_bw() +
  xlim(0,150) +
  ggtitle(label = "T. fasciculata") + 
  theme(plot.title = element_text(size = 11, face = "italic"))  +
  geom_vline(xintercept = 46.1712, linetype = "longdash", colour = "gray28") +
  scale_color_manual(values = mycolors, name = "Duplication state") +
  scale_fill_manual(values = mycolors) +
  xlab("") +
  theme(legend.title = element_text(size=9, 
                                    face="bold")) +
  theme(legend.text = element_text(size=9)) +
  theme(plot.margin = unit(c(.1,.1,.001,.1), "cm"))

p2 <- ggplot(tlei_genes, aes(x=mean_cov, colour=duplicated)) +
  geom_density() +
  theme_bw() +
  xlim(0,150) +
  ggtitle(label = "T. leiboldiana") + 
  theme(plot.title = element_text(size = 11, face = "italic"))  +
  geom_vline(xintercept = 53.34, linetype = "longdash", color = "gray28") +
  scale_color_manual(values = mycolors, name = "Duplication state") +
  scale_fill_manual(values = mycolors) +
  xlab("") + ylab("") +
  theme(plot.margin = unit(c(.1,.1,.001,.1), "cm"))

# Extract legend
legend1 <- get_legend(p1)
#Remove the legend from the box plot
p1 <- p1 + theme(legend.position="none")
p2 <- p2 + theme(legend.position="none")
pdf(paste0("FigureS1_Multi-copy_gene_assessment.pdf"), width = 12, height = 5)
plot1 <- grid.arrange(p1, p2, nrow = 1, legend1, 
                      bottom=titleGrob("Mean per-gene coverage", gp=gpar(fontsize=10)),
                      widths=c(20,20,10))
dev.off()
