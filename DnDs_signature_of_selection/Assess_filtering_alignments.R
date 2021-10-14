setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/DnDs_signature_of_selection/")
library(ggplot2)

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

ggplot(kaks_test, aes(x=filter, y=Ka.Ks, fill=alignment)) + 
         geom_boxplot() +
         ylim(c(0,5))

ggplot(kaks_test, aes(x=filter, y=Ka.Ks)) + 
  geom_boxplot() +
  ylim(c(0,5)) +
  geom_jitter(shape=16, position=position_jitter(0.2))

tapply(kaks_test$Ka.Ks, kaks_test$filter, mean)
