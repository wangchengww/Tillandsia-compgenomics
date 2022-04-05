### Gathering statistics for Structural Variant calling with SVIM ###
setwd('/Users/clara/bio-info_phd/Comparative_genomics/StrucVar/')
library(ggplot2)
library(wesanderson)

variants <- read.table('StrucVar_SVIM_raw_table.txt', header=T)
table(variants$SVTYPE)

filtered_10 <- variants[variants$QUAL >= 10,]
table(filtered_10$SVTYPE)

filtered_15 <- variants[variants$QUAL >= 15,]
table(filtered_15$SVTYPE)

ggplot(variants, aes(x=QUAL)) +
  geom_density()

### Gathering statistics for Structural Variant calling with DELLY ###
variants_5x <- read.table("StrucVar_delly5x_raw_table.txt", header = T)
table(variants_5x$SVTYPE)
variants_10x <- read.table("StrucVar_delly10x_raw_table.txt", header = T)
table(variants_10x$SVTYPE)
variants_20x <- read.table("StrucVar_delly20x_raw_table.txt", header = T)
table(variants_20x$SVTYPE)
variants_50x <- read.table("StrucVar_delly50x_raw_table.txt", header = T)
table(variants_50x$SVTYPE)

### Gathering statistics for Structural Variant calling with LUMPY ###
variants_5x_l <- read.table("StrucVar_lumpy5x_raw_table.txt", header = T)
table(variants_5x_l$SVTYPE)
variants_10x_l <- read.table("StrucVar_lumpy10x_raw_table.txt", header = T)
table(variants_10x_l$SVTYPE)
variants_20x_l <- read.table("StrucVar_lumpy20x_raw_table.txt", header = T)
table(variants_20x_l$SVTYPE)
variants_50x_l <- read.table("StrucVar_lumpy50x_raw_table.txt", header = T)
table(variants_50x_l$SVTYPE)

mycolors <- c(wes_palette("Darjeeling1")[1], wes_palette("Darjeeling1")[2], wes_palette("FantasticFox1")[3],
              wes_palette("Darjeeling1")[3], wes_palette("GrandBudapest1")[3])

summary <- read.table("illumina_SV_comparison", header = T)
summary$cov <- as.factor(summary$cov)
summary$cov <- relevel(summary$cov, "5x")
ggplot(summary, aes(x = cov, y = count, color = svtype, group = interaction(svtype, method),
                    shape = method)) + 
  geom_point(size=2) + geom_line(aes(linetype=method)) +
  scale_shape_manual(values=c(19, 1)) + 
  scale_color_manual(values = mycolors) +
  geom_hline(yintercept = 25673, color = mycolors[4], linetype = "dotdash") +
  geom_hline(yintercept = 2, color = mycolors[5], linetype = "dotdash") +
  geom_hline(yintercept = 34423, color = mycolors[2], linetype = "dotdash") +
  geom_hline(yintercept = 375, color = mycolors[3], linetype = "dotdash") +
  geom_hline(yintercept = 1372, color = mycolors[1], linetype = "dotdash") 

ggplot(summary, aes(x = cov, y = proportion, color = svtype, group = interaction(svtype, method),
                    shape = method)) + 
  geom_point(size=2) + geom_line(aes(linetype=method)) +
  scale_shape_manual(values=c(19, 1)) +
  scale_color_manual(values = mycolors) +
  geom_hline(yintercept = 41.51184413, color = mycolors[4], linetype = "dotdash") +
  geom_hline(yintercept = 0.00323389118, color = mycolors[5], linetype = "dotdash") +
  geom_hline(yintercept = 55.66011804, color = mycolors[2], linetype = "dotdash") +
  geom_hline(yintercept = 0.6063545962, color = mycolors[3], linetype = "dotdash") +
  geom_hline(yintercept = 2.218449349, color = mycolors[1], linetype = "dotdash")

overlap <- read.table("sharing_SVcallers_percentage_broad",header=T)
overlap$cov <- as.factor(overlap$cov)
overlap$cov <- relevel(overlap$cov, "5x")
ggplot(overlap, aes(x = cov, y = perc, color = sharing, group = sharing,)) + 
  geom_point(size=2) + geom_line() +
  scale_color_manual(values = mycolors)

overlap_twocallers <- read.table("sharing_SVcallers_percentage_twocallers", header=T)
overlap_twocallers$cov <- as.factor(overlap_twocallers$cov)
overlap_twocallers$cov <- relevel(overlap_twocallers$cov, "5x")
ggplot(overlap_twocallers, aes(x = cov, y = perc, color = sharing, group = sharing,)) + 
  geom_point(size=2) + geom_line() +
  scale_color_manual(values = mycolors)

agreement_threecallers <- read.table("Agreement_SVcallers_percentage_threecallers", header=T)
agreement_threecallers$cov <- as.factor(agreement_threecallers$cov)
agreement_threecallers$cov <- relevel(agreement_threecallers$cov, "5x")
ggplot(agreement_threecallers, aes(x = cov, y = perc, color = agreement, group = agreement,)) + 
  geom_point(size=2) + geom_line() +
  scale_color_manual(values = mycolors)

agreement_twocallers <- read.table("Agreement_SVcallers_percentage_twocallers", header=T)
agreement_twocallers$cov <- as.factor(agreement_twocallers$cov)
agreement_twocallers$cov <- relevel(agreement_twocallers$cov, "5x")
ggplot(agreement_twocallers, aes(x = cov, y = perc, color = agreement, group = agreement,)) + 
  geom_point(size=2) + geom_line() +
  scale_color_manual(values = mycolors)

private <- read.table("private_calls_prop", header=T)
private$cov <- as.factor(private$cov)
private$cov <- relevel(private$cov, "5x")
ggplot(private, aes(x = cov, y = perc, color = type, group = type,)) + 
  geom_point(size=2) + geom_line() +
  scale_color_manual(values = mycolors)
