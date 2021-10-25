setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/DnDs_signature_of_selection/")
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/DnDs_signature_of_selection/")
library(ggplot2)
library("tidyverse")

# Merge different alignment and filter runs into one large dataframe
kaks_strinfil <- read.delim("KaKs_stringent_filtering.txt", sep = '\t', header=T)
kaks_reldif <- read.delim("KaKs_relaxed_reldif.txt", sep = '\t', header=T)
kaks_comp <- read.delim("KaKs_relaxed_comp.txt", sep = '\t', header=T)
kaks_reldif_comp <- read.delim("KaKs_relaxed_reldif_comp.txt", sep = '\t', header=T)
kaks_ext_strinfil <- read.delim("KaKs_extensive_stringent_filtering.txt", sep = '\t', header=T)
kaks_ext_reldif <- read.delim("KaKs_extensive_relaxed_reldiff.txt", sep = '\t', header=T)
kaks_ext_comp <- read.delim("KaKs_extensive_relaxed_comp.txt", sep = '\t', header=T)
kaks_ext_reldif_comp <- read.delim("KaKs_extensive_relaxed_reldiff_comp.txt", sep = '\t', header=T)

kaks_test <- rbind(kaks_strinfil, kaks_reldif, kaks_comp, kaks_reldif_comp, kaks_ext_strinfil,
                   kaks_ext_reldif, kaks_ext_comp, kaks_ext_reldif_comp)
filter <- c(rep("stringent_filtering", 200), rep("relaxed_length_diff", 200),
            rep("relaxed_completeness", 199), rep("relaxed_length_completeness", 200),
            rep("stringent_filtering", 200), rep("relaxed_length_diff", 200),
            rep("relaxed_completeness", 200), rep("relaxed_length_completeness", 200))
alignment <- c(rep("basic", 799), rep("extensive", 800))

kaks_test <- cbind(kaks_test, filter, alignment)
kaks_test_na <- kaks_test %>% drop_na(Ka.Ks)

ggplot(kaks_test, aes(x=filter, y=Ka.Ks, fill=alignment)) + 
         geom_boxplot() +
         ylim(c(0,5))

ggplot(kaks_test_na, aes(x=filter, y=(Ka.Ks))) + 
  geom_boxplot() +
  ylim(c(0,5)) +
  geom_jitter(shape=16, position=position_jitter(0.2))

# Count NAs per filter category              
table(kaks_test[is.na(kaks_test$Ka),25]) # 10 stringent, 8 relaxed length diff, 16 relaxed completeness, 12 relaxed length and completeness
table(kaks_test[is.na(kaks_test$Ks),25]) # 10 stringent, 12 relaxed length diff, 8 relaxed completeness, 10 relaxed length and completeness

nrow(kaks_test[kaks_test$P.Value.Fisher. < 0.05 & kaks_test$filter == "stringent_filtering",]) # 216
nrow(kaks_test[kaks_test$P.Value.Fisher. < 0.05 & kaks_test$filter == "relaxed_completeness",]) # 220
nrow(kaks_test[kaks_test$P.Value.Fisher. < 0.05 & kaks_test$filter == "relaxed_length_diff",]) # 178
nrow(kaks_test[kaks_test$P.Value.Fisher. < 0.05 & kaks_test$filter == "relaxed_length_completeness",]) # 206

summary(kaks_test)
tapply(kaks_test_na$Ka.Ks, kaks_test_na$filter, mean)

checklist <- read.delim("checklist_with_diff_completeness.txt", header = T)
table(checklist$completeness)
hist(checklist$diff_per_Tfas)
ggplot(checklist, aes(x=diff_per_Tfas)) + geom_histogram(binwidth = 0.1) +
  xlim(c(-1,1))

# set up cut-off values 
breaks <- seq(from = -1, to = 1, by = 0.1)
# specify interval/bin labels
n=1
tags <- c()
for (i in 1:(length(breaks)-1)){
  tags <- c(tags, paste('[',breaks[i],']-[', breaks[i+1], ']', sep =''))
  n = n+1
}
# bucketing values into bins
group_tags <- cut(checklist$diff_per_Tfas, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
quantile(checklist$diff_per_Tfas, probs=c(.05,.15,.25,.5,.75,.85,.95))
nrow(checklist[checklist$diff_per_Tfas > -0.2 & checklist$diff_per_Tfas < 0.2 & checklist$diff_per_Tlei > -0.2 & checklist$diff_per_Tlei < 0.2,])
# 10,577
nrow(checklist[checklist$diff_per_Tfas > -0.2 & checklist$diff_per_Tfas < 0.2 & checklist$diff_per_Tlei > -0.2 & checklist$diff_per_Tlei < 0.2 & checklist$difference < 200,])
# 10,199
table(checklist[checklist$diff_per_Tfas > -0.2 & checklist$diff_per_Tfas < 0.2 & checklist$diff_per_Tlei > -0.2 & checklist$diff_per_Tlei < 0.2,7])
# 0 complete    both complete     one complete 
#      215          7241          3121
# inspect bins
summary(group_tags)
quantile(checklist$difference, c(.05,.1,.9,.95))
quantile(checklist[checklist$diff_per_Tfas > -0.2 & checklist$diff_per_Tfas < 0.2 & checklist$diff_per_Tlei > -0.2 & checklist$diff_per_Tlei < 0.2,4], c(.02,.98))
