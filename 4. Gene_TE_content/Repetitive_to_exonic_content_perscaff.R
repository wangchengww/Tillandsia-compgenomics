setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/4. Gene_TE_content/")
library(ggplot2)
library(tidyverse)
#install.packages("ggpubr")
library(ggpubr)
library(rstatix)
# Read data
tfas <- read.delim("Tillandsia_fasciculata.perscaff.genic_and_TE_content.txt", header = T, 
                   sep = "\t")
tlei <- read.delim("Tillandsia_leiboldiana.perscaff.genic_and_TE_content.txt", header = T, 
                   sep = "\t")

tfas$repetitive_content <- (tfas$repetitive_content)*0.01
tfas$repetitive_bases <- tfas$total_length * tfas$repetitive_content
tfas$ratio <- tfas$total_exonic_length / tfas$repetitive_bases
tfas$ratio2 <- tfas$repetitive_bases / tfas$total_exonic_length

tlei$repetitive_content <- (tlei$repetitive_content)*0.01
tlei$repetitive_bases <- tlei$total_length * tlei$repetitive_content
tlei$ratio <- tlei$total_exonic_length / tlei$repetitive_bases
tlei$ratio2 <- tlei$repetitive_bases / tlei$total_exonic_length

ratios <- as.data.frame(cbind(c(tfas$Scaffold, tlei$Scaffold), 
                              c(tfas$ratio2, tlei$ratio2),
                              c(rep("T.fasciculata", 25), rep("T.leiboldiana", 26))))
ratios$rep_to_exon_ratio <- as.numeric(ratios$rep_to_exon_ratio)
colnames(ratios) <- c("scaff", "rep_to_exon_ratio", "species")

ggplot(ratios, aes(x = species,y = rep_to_exon_ratio)) + geom_boxplot()
ggplot(ratios, aes(x = rep_to_exon_ratio, color = species)) + geom_histogram(binwidth = 1)

# Test for normality in each group
ratios %>%
  group_by(species) %>%
  shapiro_test(rep_to_exon_ratio)
# A tibble: 2 × 4
# species       variable          statistic     p
# <chr>         <chr>                 <dbl> <dbl>
# 1 T.fasciculata rep_to_exon_ratio     0.974 0.745 #Is normal
# 2 T.leiboldiana rep_to_exon_ratio     0.953 0.278 #Is normal

# Draw a qq plot by group
ggqqplot(ratios, x = "rep_to_exon_ratio", facet.by = "species")
# The ratios are normally disributed in both species

# Levene test, equality of variances. If equal, we can use T-test, otherwise Welch's test
ratios %>% levene_test(rep_to_exon_ratio ~ species)
# A tibble: 1 × 4
# df1   df2 statistic           p
# <int> <int>     <dbl>       <dbl>
#  1     1    49      32.5 0.000000678 #Variances are not equal --> Welch's test

# Welch's test, is the repetitive-to-exonic content ratio per scaffold significantly
# different between Tfas and Tlei? (((Yes)))
stat.test <- ratios %>% 
  t_test(rep_to_exon_ratio ~ species) %>%
  add_significance()
stat.test
# A tibble: 1 × 9
#.y.               group1        group2           n1    n2 statistic    df        p p.signif
#* <chr>             <chr>         <chr>         <int> <int>     <dbl> <dbl>    <dbl> <chr>   
#  1 rep_to_exon_ratio T.fasciculata T.leiboldiana    25    26     -4.00  27.6 0.000433 *** 
mean(ratios$rep_to_exon_ratio)
mean(tfas$ratio2) #27.6
mean(tlei$ratio2) #46.8
mean(tfas$total_exonic_length)
mean(tlei$total_exonic_length)
