### Assessment of gene duplications in Tfas ###

### Here I check whether inferred gene duplication in Tillandsia fasciculata are legitimate or not
### based on median coverage calculations after mapping raw illumina data, 50x, to the reference of 
### Tfas. The idea is that genuine duplications will maintain stable coverage which is similar to
### the background, whereas false duplications will have lower coverage. 

library(ggplot2)
library(plyr)
library(dplyr)
library(wesanderson)
wes_palette("Darjeeling1")[1]
# Set colours
mycolors <- c(wes_palette("Darjeeling1")[1], wes_palette("Darjeeling1")[2], wes_palette("FantasticFox1")[3],
              wes_palette("Darjeeling1")[3], wes_palette("Rushmore1")[3], wes_palette("GrandBudapest1")[3])

setwd("/Users/clara/bio-info_phd/Comparative_genomics/Gene_family_evolution/")
tfas_genes <- read.table("Tfas_pergene_mediancov_and_orthoinfo.txt", header = T)
mean(tfas_genes$median_cov) # 47.3
mean(tfas_genes$mean_cov) # 47.9

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
mu <- ddply(tfas_genes, "dup", summarise, grp.mean=mean(mean_cov))
head(mu)

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


### Establish a cutoff to separate "real" from "faulty" gene models
## Using the package "cutoff", we can establish a cut-off value that separates the bimodal
## distribution into two unimodal, normal distributions. The cut-off value represents the point
## where the likelihood of a datapoint belonging only to one of the two distributions is high 
## enough that we can call it with confidence. This will be a slightly conservative value, but
## allows us to separate the two groups based on statistics and probability. 

devtools::install_github("choisy/cutoff")
install.packages("cli")
install.packages("bbmle")
library(cutoff)
??cutoff
tfas_multicopy_mean_cov <- as.vector(tfas_genes[tfas_genes$duplicated == 'multi-copy', 3])
tfas_multicopy_mean_cov_df <- tfas_genes[tfas_genes$duplicated == 'multi-copy',]
tfas_multicopy_mean_cov_subset <- tfas_multicopy_mean_cov[tfas_multicopy_mean_cov < 100]
tfas_multicopy_FMM <- em(tfas_multicopy_mean_cov_subset, 'normal', 'normal')
hist(tfas_multicopy_mean_cov_subset,10000,F,xlim=c(0,100),xlab="median coverage",ylab="density",main=NULL,col="grey")
lines(tfas_multicopy_FMM,lwd=1.5,col="red")

cut_off <- cutoff(tfas_multicopy_FMM)
cut_off # Genes with median cov < 34.5 belong to "faulty" category

faulty_genes <- tfas_multicopy_mean_cov_df[tfas_multicopy_mean_cov_df$mean_cov < 34.5,] # 4881 genes
table(faulty_genes$Tfas_count)
orthogroups_w_faulty_genes <- tfas_genes %>%
  filter(og_id %in% faulty_genes$og_id)

high_cov_genes <- tfas_genes[tfas_genes$median_cov > 100,]
write.table(high_cov_genes, file = "100x_genes_Tfas.txt")
write.table(faulty_genes, file = "faulty_gene_families_35x.txt")
write.table(orthogroups_w_faulty_genes, file = "ogs_with_faulty_genes_35x.txt")

# Correcting the gene familz sizes based on coverage differences
corr_Tfas <- data.frame()
for (i in unique(orthogroups_w_faulty_genes$og_id)){
  orthogroup <- orthogroups_w_faulty_genes[orthogroups_w_faulty_genes$og_id == i,]
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
rownames(corr_Tfas) <- c(1:1981)
colnames(corr_Tfas) <- c("og_id", "correction_factor", "corr_Tfas_count","old_Tfas_count", "total_mean_cov", "expected_mean_cov")

write.table(corr_Tfas, file = "corrected_family_sizes_Tfas.txt")
###########################################################################################
### Same for Tleiboldiana

tlei_genes <- read.table("Tlei_pergene_mediancov_and_orthoinfo.txt", header = T)
mean(tlei_genes$median_cov) # 66.9
mean(tlei_genes$mean_cov) # 66.6

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
  geom_vline(xintercept = 46.1712, linetype = "longdash", colour = "gray28") +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)

