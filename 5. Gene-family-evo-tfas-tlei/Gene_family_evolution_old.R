## Studying Gene Family evolution between T. fasciculata and T. leiboldiana using GLMs

## Goal:
# We have inferred orthogroups (sets of genes between species which descend from a
# common ancestor) between T.fasciculata and T. leiboldiana. The number of genes per
# species in each orthogroup varies. Most orthogroups contain one gene of both species,
# but we have a number of orthogroups where the numbers vary. By fitting a model to
# the distribution of counts (linear regression or GLM), we may be able to
# identify orthogroups with significant differences in gene count between both species.
# This would be a first, statistically demonstrated bit of evidence for gene family evolution
# between both species. One could follow up with GO term enrichment in these outlier
# orthogroups to find out the functions of these evolving gene families.

# 1. Visualizing and exploring the data

## setwd and import data
setwd("/Users/clara/bio-info_phd/Comparative_genomics/Gene_family_evolution//")
genes <- read.delim("orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt", sep = "\t", header = F)
counts <- read.table("orthogroup_counts.txt", sep = '\t')
colnames(counts) <- c("og_id", "Acom", "Tfas", "Tlei")
summary(counts)
# Here we select orthogroups which contain at least one gene in both species. I do this
# mostly to avoid "nonsense" orthogroups, as species-specific annotations are often TEs
# or didn't blast during annotation, so we can't safely say they reflect real genes.
counts_Tfas_Tlei <- counts[counts$Tfas != 0 & counts$Tlei != 0,]
counts_Tfas_Tlei$diff <- counts_Tfas_Tlei$Tfas - counts_Tfas_Tlei$Tlei

d <- counts_Tfas_Tlei[counts_Tfas_Tlei$diff > 8 | counts_Tfas_Tlei$diff < -8,]
write.table(d, file = "orthogroups_size_difference_8.txt", sep = "\t")
# I remove here an additional orthogroup which is a strong outlier in Tlei and has no
# bkast results from the annotation.
counts_Tfas_Tlei <- counts_Tfas_Tlei[counts_Tfas_Tlei$og_id != "OG0000015",]
counts_Tfas_Tlei_multi <- counts_Tfas_Tlei[!(counts_Tfas_Tlei$Tfas == 1 & counts_Tfas_Tlei$Tlei == 1),]
onetoone <-  counts[counts$Tfas == 1 & counts$Tlei == 1,]

#Check the skewness of our data
library(e1071)
skewness(counts_Tfas_Tlei$Tfas) #23.34
skewness(counts_Tfas_Tlei$Tlei) #11.7
# Check skewness after removing all 1-1 orthogroups
skewness(counts_Tfas_Tlei_multi$Tfas) #6.59
skewness(counts_Tfas_Tlei_multi$Tlei) #5.6
# Check skewness after a log transformation
skewness(log(counts_Tfas_Tlei$Tfas)) #3.54
skewness(log(counts_Tfas_Tlei$Tlei)) #4.02

# Visualizing distribution of the counts
library(reshape2)
library(ggplot2)
counts_Tfas_Tlei_t <-melt(counts_Tfas_Tlei[,c(1,3,4)], id.vars = c("og_id"))
colnames(counts_Tfas_Tlei_t) <- c('og_id', 'species', 'gene_count')
ggplot(counts_Tfas_Tlei_t, aes(gene_count, fill = species)) +
  geom_histogram(binwidth = 1) + facet_grid(species ~ ., margins = TRUE, scales = "free")

#Make a scatterplot and fit a linear regression line
library(ggplot2)
ggplot(counts_Tfas_Tlei_multi, aes(x=Tfas, y=Tlei)) + geom_point(size = .5) +
  #geom_smooth(method=lm, se=F) +
  geom_abline(intercept = 0, slope = 1) +
  ylim(0,85) +
  xlim(0,85) +
  labs(title = "Per-species gene counts in multi-copy orthogroups") +
  xlab(label = "T. leiboldiana") +
  ylab(label = "T. fasciculata")
  theme_bw()
fit<- lm(Tlei ~ Tfas, data=counts_Tfas_Tlei)
summary(fit)
plot(fit)

# As is clear from the skewdness and the data visualization, the counts are not
# normally distributed. Even though the linear regression is strongly significant,
# the R2 value is very low (0.089) and we are violating the conditions for linear
# regression. This is very clear when plotting the diagnostic plots of the model
# plot(fit).

# Interestingly, the plot shows us an assymetric relationship of gene counts between
# the two species. It seems that generally, T. fasciculata tends to harbour much
# larger numbers of genes in orthogroups than T. leiboldiana.

# We can try to normalize the data by performing a log transformation, as the
# skewness is reduced when logtransforming the data (but still skewed, as the value
# is > 1).
ggplot(data = counts_Tfas_Tlei, aes(x = Tfas, y = Tlei)) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method=lm, se=T) +
  geom_abline(intercept = 0, slope = 1)
fit_log.model <- lm(log1p(Tlei) ~ log1p(Tfas), data = counts_Tfas_Tlei)
summary(fit_log.model)
plot(fit_log.model)

# Though the R2 is higher in this model, we still see very abnormal diagnostic plots,
# showing that logtransforming our data is not helping reaching the assumptions
# necessary for a linear regression.

# Clearly, transforming the data doesn't correct for the skewdness. Linear regression
# is therefore not an option here. Generalized Linear Models can be used for data
# with a distribution different from normal.

# But, what is the distribution of our data? Our data is of the "count" type, we
# have the count of genes for each species in each orthogroup. There are no
# negative values and only integers. For this kind of data, the poisson or negative
# binomial distributions can apply.

# The difference between these two is the relationship of the mean to the variance.
# In a poisson distribution, the mean = variance, whereas in negative binomial, the
# mean < variance. In the latter, the data is very dispersed.

# 2. Fitting generalized linear models to our data

#Let's see how the mean relates to the variance in our data:

mean(counts_Tfas_Tlei$Tfas)   # 1.36
var(counts_Tfas_Tlei$Tfas)    # 4.56
mean(counts_Tfas_Tlei$Tlei)   # 1.17
var(counts_Tfas_Tlei$Tlei)    # 0.63

intercept <- mean(counts_Tfas_Tlei$Tfas)
# It seems that for T.fasciculata, the variance is indeed larger than the mean,
# but for T. leiboldiana this doesn't seem to be the case. Therefore, I computed
# both a poisson GLM and a negative binomial (NB) one and compared them with the
# likelihood-ratio test:

library(MASS)
library(lmtest)
fit.nb <- glm.nb(gene_count ~ species, data = counts_Tfas_Tlei_t)
fit.poisson <- glm(gene_count ~ species, family = "poisson", data = counts_Tfas_Tlei_t)
summary(fit.poisson)
summary(fit.nb)
# test if there's a significant difference if removing "species"
m2 <- update(fit.nb, . ~ . - species)
summary(m2)
a <- anova(fit.nb, m2)
# Compare poisson to negbin
lrtest(fit.poisson, fit.nb)
pchisq(fit.nb$deviance, df=fit.nb$df.residual, lower.tail=FALSE)
pchisq(fit.poisson$deviance, df=fit.poisson$df.residual, lower.tail=FALSE)

ggplot(data = counts_Tfas_Tlei, aes(x = Tfas, y = Tlei)) +
  geom_point() +
  geom_smooth(method=glm.nb, se=T)

## Alternatively, one can formally test whether there is dispersion in the model:
library(AER)
deviance(fit.poisson)/fit.poisson$df.residual
dispersiontest(fit.poisson)
# The p-value is < 0.05, meaning the alt hypothesis is true. So we have overdispersion.
lrtest()
# Likelihood ratio test
# Model 1: gene_count ~ species
# Model 2: gene_count ~ species
# Df LogLik Df  Chisq Pr(>Chisq)
# 1   2 -42998
# 2   3 -42614  1 768.89  < 2.2e-16 ***
#
# So, the likelihood ratio test argues that our negative binomial model fits better
# than the model based on a Poisson distribution. The respective AIC of both models
# supports this (85234 in nb and 86001 in poisson). The deviance residuals are also
# smaller in nb than in poisson.

# We can have a closer look at the distribution for with the countreg package, and
# simulate data with the models parameters to see how well it resembles to the
# observed data:

#install.packages("countreg", repos="http://R-Forge.R-project.org")
library(countreg)
countreg::rootogram(fit.nb)
countreg::rootogram(fit.poisson)

library(magrittr)
op <- par(mfrow=c(1,2))
set.seed(1)
counts_Tfas_Tlei_t$gene_count %>% `[`(counts_Tfas_Tlei_t$species=="Tfas") %>%
  table() %>% barplot(main = "Observed Tfas")
rnbinom(n = 16424, size = fit.nb$theta, mu = exp(coef(fit.nb)[1])) %>%
  table() %>%  barplot(main = "Simulated Tfas")

counts_Tfas_Tlei_t$gene_count %>% `[`(counts_Tfas_Tlei_t$species=="Tlei") %>%
  table() %>% barplot(main = "Observed Tlei")
rnbinom(n = 16424, size = fit.nb$theta, mu = exp(sum(coef(fit.nb)))) %>%
  table() %>%  barplot(main = "Simulated Tlei")
par(op)

## Plot residuals
res <- (residuals(fit.nb, type="deviance"))
hist(res, breaks = 50)
fitted <- as.data.frame(predict(fit.nb, type = "response"))
plot((predict(fit.nb, type = "response")), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)
?predict()
library(AER)
deviance(fit.poisson)/fit.poisson$df.residual
dispersiontest(fit.poisson)

# After discussing a bit with Aglaia and Marta, I decided to delve a bit into the zero-inflated
# models. According to Marta, there should be a way of modifying existing models towards one-inflated.
# I first try a few R packages that already seem to have this incorporated:
require(pscl)
require(MASS)
require(boot)
# I subtract each count with 1 to obtain zero-inflation instead of one inflation.
# This was just meant as a test to see if zero-inflation actually does improve anything.
counts_Tfas_Tlei_t$gene_count_zero <- counts_Tfas_Tlei_t$gene_count - 1
fit.zerinfl <- zeroinfl(gene_count_zero ~ species, data = counts_Tfas_Tlei_t,
                              dist = "negbin")
m2 <- glm.nb(gene_count_zero ~ species, data = counts_Tfas_Tlei_t)
summary(fit.zerinfl)
lrtest(m2, fit.zerinfl)
vuong(fit.zerinfl, m2)

# Based on the above tests, the zero-inflated model in fact does not improve the fit. So this
# doesn't seem the right way to go. And in any case, modifying the data like this is a very
# sketchy thing to do

# I also tested out VGAM's one-inflated poisson model, but this model doesn't have a deviance function.
# Therefore, I can't plot deviance residuals against the predictors.
library(VGAM)
library(VGAMdata)
fit.oneinfl.pois <- vglm(gene_count ~ species, oipospoisson, data = counts_Tfas_Tlei_t)
summary(fit.oneinfl.pois)
# I also computed a normal negative binomial model (until now the best fit) in VGAM for comparison
fit.vglm.nb <- vglm(gene_count ~ species, negbinomial(), data = counts_Tfas_Tlei_t)
summary(fit.vglm.nb)
lrtest_vglm(fit.vglm.nb, fit.oneinfl.pois)

# Model 1: gene_count ~ species
# Model 2: gene_count ~ species
# Df LogLik Df Chisq Pr(>Chisq)
# 1 65693 -42614
# 2 65692 -21311 -1 42606  < 2.2e-16 ***
# This seems to indicate that our one-inflated poisson model better fits our data than the classical
# negative binomial.
plotvgam(fit.oneinfl.pois, residuals = T, type.residuals = "working", se = T)
plot(predict(fit.oneinfl.pois), residuals(fit.oneinfl.pois, type="pearson"))
pchisq(fit.oneinfl.pois, df=fit.poisson$df.residual)

## Trying the same thing but without the 1-1 orthologs
#
library(reshape2)
var(counts_Tfas_Tlei_multi$Tfas)
var(counts_Tfas_Tlei_multi$Tlei)
ggplot(counts_Tfas_Tlei_multi, aes(x=Tfas, y=Tlei)) + geom_point(size = .5) +
  geom_smooth(method=lm, se=T) +
  geom_abline(intercept = 0, slope = 1)

counts_Tfas_Tlei_multi_t <-melt(counts_Tfas_Tlei_multi[,c(1,3,4)], id.vars = c("og_id"))
colnames(counts_Tfas_Tlei_multi_t) <- c('og_id', 'species', 'gene_count')
ggplot(counts_Tfas_Tlei_multi_t, aes(gene_count, fill = species)) +
  geom_histogram(binwidth = 1) + facet_grid(species ~ ., margins = TRUE, scales = "free")
with(counts_Tfas_Tlei_multi_t, tapply(gene_count, species, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), var(x))
}))

fit.nb_multi <- glm.nb(gene_count ~ species, data = counts_Tfas_Tlei_multi_t)
fit.poisson_multi <- glm(gene_count ~ species, family = "poisson", data = counts_Tfas_Tlei_multi_t)
summary(fit.nb_multi)
summary(fit.poisson_multi)
lrtest(fit.poisson_multi, fit.nb_multi)

countreg::rootogram(fit.nb_multi)
countreg::rootogram(fit.poisson_multi)

## Plot residuals
res <- as.data.frame(residuals(fit.nb_multi, type="deviance"))
hist(res$`residuals(fit.nb_multi, type = "deviance")`)
plot((predict(fit.nb_multi)), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)
plot(density(resid(m0, type='deviance')))
plot(resid(fit.nb, type = "pearson") ~ fitted(fit.nb))

## Trying the same thing but with all counts (also species specific orthogroups)
library(reshape2)
mean(counts$Tfas)   # 1.36
var(counts$Tfas)    # 32.33
mean(counts$Tlei)   # 1.17
var(counts$Tlei)
counts_t <-melt(counts[,c(1,3,4)], id.vars = c("og_id"))
colnames(counts_t) <- c('og_id', 'species', 'gene_count')
fit.nb_all <- glm.nb(gene_count ~ species, data = counts_t)
fit.poisson_all <- glm(gene_count ~ species, family = "poisson", data = counts_t)
summary(fit.nb_all)
summary(fit.poisson_all)
lrtest(fit.poisson_all, fit.nb_all)
countreg::rootogram(fit.nb_all)
countreg::rootogram(fit.poisson_all)

library(magrittr)
op <- par(mfrow=c(1,2))
set.seed(1)
counts_t$gene_count %>% `[`(counts_t$species=="Tfas") %>%
  table() %>% barplot(main = "Observed Tfas")
rnbinom(n = 19101, size = fit.nb_all$theta, mu = exp(coef(fit.nb_all)[1])) %>%
  table() %>%  barplot(main = "Simulated Tfas")

counts_t$gene_count %>% `[`(counts_t$species=="Tlei") %>%
  table() %>% barplot(main = "Observed Tlei")
rnbinom(n = 19101, size = fit.nb_all$theta, mu = exp(sum(coef(fit.nb_all)))) %>%
  table() %>%  barplot(main = "Simulated Tlei")
par(op)
