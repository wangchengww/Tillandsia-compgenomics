setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/3b. Chromosome 14 Tlei/')

pacman::p_load("ggplot2", "gtools")


data <- read.table("Tlei_exonic-to-repetitive_content_perscaff.txt", sep = "\t", header = T)
data$Scaffold <- factor(data$Scaffold, levels=mixedsort(unique(data$Scaffold)))

ggplot(data, aes(x=Scaffold, y=total_exonic_length)) + geom_point()
ggplot(data, aes(x=Scaffold, y=proportion_exonic)) + geom_point()
ggplot(data, aes(x=Scaffold, y=ratio2)) + geom_point()

degenes <- read.delim("orthogroup_info_DEGENES_orthologs_Tlei.txt",
                      header = F, sep = "\t")
degenes <- degenes[startsWith(degenes$V1,"Tlei"),]
degenes <- degenes[,c(1,2)]
freq <- as.data.frame(table(degenes$V2))
freq$Var1 <- factor(freq$Var1, levels=mixedsort(unique(freq$Var1)))
colnames(freq) <- c("Scaffold", "DE_count")
data <- merge(data, freq, by="Scaffold")

allgenes <- read.delim("orthogroup_info_ALLGENES_Tlei.txt", header = F, sep = "\t")
allgenes <- allgenes[startsWith(allgenes$V1,"Tlei"),]
allgenes <- allgenes[,c(1,2)]
freq2 <- as.data.frame(table(allgenes$V2))
colnames(freq2) <- c("Scaffold", "all_genecount")
data <- merge(data, freq2, by="Scaffold")
data$proportion_degenes <- data$DE_count/data$all_genecount

ggplot(data, aes(x=Scaffold, y=DE_count)) + geom_point()
