library(ggplot2)
setwd("/Users/clara/Documents/GitHub/Tillandsia-compgenomics/6. DnDs_signature_of_selection/")

dnds_all <- read.delim("pairwise_dnds.codeml.allgenes.output", sep = "\t", header = T)
dnds_possel_sign <- dnds_all[dnds_all$dN.dS > 1 & dnds_all$pvalue < 0.05,]
dnds_possel <- dnds_all[dnds_all$dN.dS > 1,]
dnds_sign <- dnds_all[dnds_all$pvalue < 0.05,]

c <- c()
for (i in 1:nrow(dnds_all)){
  if (dnds_all[i,8] <= 0.05){
    c <- c(c, "sign")
  } else {
    c <- c(c, "non")
  }
}
dnds_all$sign <- c

ggplot(dnds_sign) + aes(x = dN.dS) + geom_density() + xlim(c(-1,3))
ggplot(dnds_all) + aes(x = dN.dS, color = sign) + geom_density() + xlim(c(-1,3))

mean(dnds_sign$dN.dS)
mean(dnds_sign$dN)
mean(dnds_sign$dS)
