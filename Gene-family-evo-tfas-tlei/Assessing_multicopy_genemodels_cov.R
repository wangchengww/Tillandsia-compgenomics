### Assessment of gene duplications in Tfas ###

### Here I check whether inferred gene duplication in Tillandsia fasciculata are legitimate or not
### based on median coverage calculations after mapping raw illumina data, 50x, to the reference of
### Tfas. The idea is that genuine duplications will maintain stable coverage which is similar to
### the background, whereas false duplications will have lower coverage.

library(ggplot2)
library(plyr)
library(dplyr)
library(wesanderson)

# Set colours categories
mycolors <- c(wes_palette("Darjeeling1")[1], wes_palette("Darjeeling1")[2], wes_palette("FantasticFox1")[3],
              wes_palette("Darjeeling1")[3], wes_palette("Rushmore1")[3], wes_palette("GrandBudapest1")[3])

setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
tfas_genes <- read.table("Tfas_pergene_mediancov_and_orthoinfo.txt", header = T)

# Removing all genes from chloroplastic, ribosomal or mitochondrial genes
mito_og <- read.table("mito_plastid_ribo_OGs.txt", header = F, sep = "\t")
colnames(mito_og) <- c("og_id")
tfas_no_mito <- tfas_genes %>%
  filter(!(og_id %in% mito_og$og_id))

mean(tfas_genes$median_cov) # 47.5
mean(tfas_genes$mean_cov) # 48.04

# Divide genes into categories based on copy number
dup <- c()
for (i in 1:nrow(tfas_genes)){
  if (tfas_genes[i,6] == 1){
    dup <- c(dup, "single-copy")
  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,5] == 0 && tfas_genes[i,7] == 0){
    dup <- c(dup, "unique-multi-copy")
  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,6] > tfas_genes[i,5] && tfas_genes[i,6] > tfas_genes[i,7]){
    dup <- c(dup, "largest-multi-copy")
  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,6] < tfas_genes[i,5] && tfas_genes[i,6] < tfas_genes[i,7]){
    dup <- c(dup, "smallest-multi-copy")
  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,6] == tfas_genes[i,5] && tfas_genes[i,6] == tfas_genes[i,7]){
    dup <- c(dup, "ancestral-multi-copy")
  } else {
    dup <- c(dup, "middle-multi-copy")
  }
}
tfas_genes$duplicated <- dup

# simplified categories
dup <- c()
for (i in 1:nrow(tfas_genes)){
  if (tfas_genes[i,6] == 1 && tfas_genes[i,5] == 1 && tfas_genes[i,7] ==1){
    dup <- c(dup, "ancestral-single-copy")
  } else if (tfas_genes[i,6] == 1){
    dup <- c(dup, "single-copy")
  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,5] == 0 && tfas_genes[i,7] == 0){
    dup <- c(dup, "unique-multi-copy")
  } else if (tfas_genes[i,6] > 1 && tfas_genes[i,6] == tfas_genes[i,5] && tfas_genes[i,6] == tfas_genes[i,7]){
    dup <- c(dup, "ancestral-multi-copy")
  } else {
    dup <- c(dup, "multi-copy")
  }
}
tfas_genes$duplicated <- dup

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

sd(tfas_no_mito[tfas_no_mito$duplicated == "ancestral-single-copy", 3])
mean(tfas_no_mito[tfas_no_mito$duplicated == "ancestral-single-copy", 3])
quantile(tfas_no_mito[tfas_no_mito$duplicated == "ancestral-single-copy", 3], c(.01, .05, .1, .25, .5,.75, .9, .95, .99)) 
#1%       5%      10%      25%      50%      75%      90%      95%      99% 
#15.10657 22.41436 27.21449 34.50790 39.82849 43.67525 46.54118 48.44374 69.23726

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

ggplot(tfas_no_mito, aes(x=mean_cov, color=duplicated)) +
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
high_cov_genes <- tfas_genes[tfas_genes$mean_cov > 100,]
write.table(high_cov_genes, file = "100x_genes_Tfas.txt")

# Recording all genes with a mean coverage < 34.5 to correct gene family sizes
#write.table(faulty_genes, file = "faulty_gene_families_35x.txt")
-----------------
  
# Correcting the gene family sizes based on coverage differences
# Select all orthogroups containing faulty genes
# library(dplyr)
# orthogroups_w_faulty_genes <- tfas_genes %>%
#  filter(og_id %in% faulty_genes$og_id)
# This resulted in 6861 genes and 1981 orthogroups

# Now I select ALL multicopy genes - this immediately captures the full orthogroup as well
tfas_multicopy_genes <- tfas_genes[tfas_genes$duplicated == 'multi-copy',] # 8400 genes, 2632 orthogroups
tfas_no_mito_multicopy_genes <- tfas_no_mito[tfas_no_mito$duplicated == 'multi-copy',] # 7348 genes, 2353 orthogroups

lower_thresh=34.5079# 25th quantile of the average coverage of ancestral single copy genes to accont for variation
upper_thresh=43.6752 # 75th quantile of the average coverage of ancestral single copy genes
corr_Tfas_no_mito <- data.frame()
for (i in unique(tfas_no_mito_multicopy_genes$og_id)){
  orthogroup <- tfas_no_mito_multicopy_genes[tfas_no_mito_multicopy_genes$og_id == i,]
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
  corr_Tfas_no_mito <- rbind(corr_Tfas_no_mito, corr_family_size)
}
colnames(corr_Tfas_no_mito) <- c("og_id", "correction_factor", "corr_Tfas_count","old_Tfas_count", "total_mean_cov", "expected_mean_cov", "Lower_thresh", "Upper_thresh")

corr_Tfas <- data.frame()
for (i in unique(tfas_multicopy_genes$og_id)){
  orthogroup <- tfas_multicopy_genes[tfas_multicopy_genes$og_id == i,]
  nr_genes <- as.integer(orthogroup[1,6])
  total_mean_cov <- sum(orthogroup$mean_cov)
  expected_mean_cov <- nr_genes * 46.1712
  correction_factor <- total_mean_cov/expected_mean_cov
  new_size = round((nr_genes * correction_factor), digits = 0)
  if (correction_factor > 1){
    new_size = nr_genes
  }
  if (total_mean_cov < 46.1712){
    new_size = 1
  }
  corr_family_size <- cbind(i, as.numeric(correction_factor), as.integer(new_size),
                            as.integer(nr_genes), as.numeric(total_mean_cov), expected_mean_cov)
  corr_Tfas <- rbind(corr_Tfas, corr_family_size)
}
colnames(corr_Tfas) <- c("og_id", "correction_factor", "corr_Tfas_count","old_Tfas_count",
                         "total_mean_cov", "expected_mean_cov")


increase <- corr_Tfas[as.integer(corr_Tfas$corr_Tfas_count) > as.integer(corr_Tfas$old_Tfas_count),]
increase_plastid <- increase %>%
  filter(!(og_id %in% mito_og$og_id))


write.table(corr_Tfas, file = "corrected_family_sizes_Tfas.txt")

###########################################################################################
### Same for Tleiboldiana
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
tlei_genes <- read.table("Tlei_pergene_mediancov_and_orthoinfo.txt", header = T)
mean(tlei_genes$median_cov) # 67.04
mean(tlei_genes$mean_cov) # 66.73

dup <- c()
for (i in 1:nrow(tlei_genes)){
  if (tlei_genes[i,7] == 1){
    dup <- c(dup, "single-copy")
  } else if (tlei_genes[i,7] > 1 && tlei_genes[i,5] == 0 && tlei_genes[i,6] == 0){
    dup <- c(dup, "unique-multi-copy")
  } else if (tlei_genes[i,7] > 1 && tlei_genes[i,7] > tlei_genes[i,5] && tlei_genes[i,7] > tlei_genes[i,6]){
    dup <- c(dup, "largest-multi-copy")
  } else if (tlei_genes[i,7] > 1 && tlei_genes[i,7] < tlei_genes[i,5] && tlei_genes[i,7] < tlei_genes[i,6]){
    dup <- c(dup, "smallest-multi-copy")
  } else if (tlei_genes[i,7] > 1 && tlei_genes[i,7] == tlei_genes[i,5] && tlei_genes[i,7] == tlei_genes[i,6]){
    dup <- c(dup, "ancestral-multi-copy")
  } else {
    dup <- c(dup, "middle-multi-copy")
  }
}
tlei_genes$duplicated <- dup

# simplified
dup <- c()
for (i in 1:nrow(tlei_genes)){
  if (tlei_genes[i,7] == 1 && tlei_genes[i,5] == 1 && tlei_genes[i,6] ==1){
    dup <- c(dup, "ancestral-single-copy")
  } else if (tlei_genes[i,7] == 1){
    dup <- c(dup, "single-copy")
  } else if (tlei_genes[i,7] > 1 && tlei_genes[i,5] == 0 && tlei_genes[i,6] == 0){
    dup <- c(dup, "unique-multi-copy")
  } else if (tlei_genes[i,7] > 1 && tlei_genes[i,7] == tlei_genes[i,5] && tlei_genes[i,7] == tlei_genes[i,6]){
    dup <- c(dup, "ancestral-multi-copy")
  } else {
    dup <- c(dup, "multi-copy")
  }
}
tlei_genes$duplicated <- dup

mu <- ddply(tlei_genes, "dup", summarise, grp.mean=mean(mean_cov))
head(mu)

ggplot(tlei_genes, aes(x=mean_cov, colour=duplicated)) +
  geom_density() +
  xlim(0,150) +
  geom_vline(xintercept = 53.34, linetype = "longdash", color = "gray28") +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)

ggplot(tlei_genes, aes(x=median_cov, color=dup)) +
  geom_density() +
  xlim(0,150) +
  geom_vline(xintercept = 53.34, linetype = "longdash", colour = "gray28") +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)

# Because the distribution here is unimodal, with just a small shoulder of faulty genes at low coverage,
# we can't separate the shoulder in the same way as we did for Tfas. So, I decided to run corrections
# on all multicopy genes.
# Correcting the gene family sizes based on coverage differences for all multicopy genes
Tlei_multicopy <- tlei_genes[tlei_genes$duplicated == "multi-copy",]
corr_Tlei <- data.frame()
for (i in unique(Tlei_multicopy$og_id)){
  orthogroup <- Tlei_multicopy[Tlei_multicopy$og_id == i,]
  nr_genes <- as.integer(orthogroup[1,7])
  total_mean_cov <- sum(orthogroup$mean_cov)
  expected_mean_cov <- nr_genes * 53.34
  correction_factor <- total_mean_cov/expected_mean_cov
  new_size = round((nr_genes * correction_factor), digits = 0)
  if (correction_factor > 1){
    new_size = nr_genes
  }
  if (total_mean_cov < 53.34){
    new_size = 1
  }
  corr_family_size <- cbind(i, as.numeric(correction_factor), as.integer(new_size),
                            as.integer(nr_genes), as.numeric(total_mean_cov), expected_mean_cov)
  corr_Tlei <- rbind(corr_Tlei, corr_family_size)
}
rownames(corr_Tlei) <- c(1:1450)
colnames(corr_Tlei) <- c("og_id", "correction_factor", "corr_Tlei_count","old_Tlei_count", "total_mean_cov", "expected_mean_cov")

write.table(corr_Tlei, file = "corrected_family_sizes_Tlei.txt")
