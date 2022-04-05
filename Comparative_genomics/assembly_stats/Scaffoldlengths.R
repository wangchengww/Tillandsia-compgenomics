## Script to display Scaffolds lengths of T.lei and T.fas ##

# Set working directory, load libraries
setwd("/Users/clara/bio-info_phd/assembly_stats/")
library(ggplot2)
library(wesanderson)

# Read raw files
lei <- read.table("scaffoldLengths_Tlei_new_assembly.txt") 
fas <- read.table("length_fas")

# Select 25 largest scaffolds
lei <-sort(lei$V2, decreasing = T)
fas$V1 <- sort(fas$V1, decreasing = T)
lei <- cbind(c(1:10433), lei)
colnames(lei) <- c("scaffold", "length_lei") 
fas <- cbind(c(1:2321), fas)
colnames(fas) <- c("scaffold", "length_fas")
lei_25 <- as.data.frame(lei[1:25,])
fas_25 <- fas[1:25,]

# Plot individual line graphs of decreasing scaffold size
ggplot(lei_25, aes(x = scaffold, y = length_lei)) + geom_line()
ggplot(fas_25, aes(x = scaffold, y = length_fas)) + geom_line()

# Create stacked dataset for histogram and density plot
both_25 <- as.data.frame(cbind(lei_25$length_lei, fas_25$length_fas))
colnames(both_25) <- c("T.leiboldiana", "T.fasciculata")
both_25 <- stack(both_25)
both_25 <- cbind(c(1:25), both_25)
colnames(both_25) <- c("scaffold", "length", "species")

# Create density plot of size distribution of 25 largest scaffolds
ggplot(both_25, aes(x = length, colour = species)) + 
  geom_density(aes(y=..scaled..)) +
  scale_x_continuous(name="Scaffold length", breaks=seq(0, 100000000, 10000000))

# Create histogram of size distribution of 25 largest scaffolds
ggplot(both_25, aes(x = length, fill = species)) + 
  geom_histogram(position="dodge", binwidth = 5000000) +
  scale_x_continuous(name="Scaffold size (Mb)", 
                     breaks=seq(0, 70000000, 10000000),
                     labels = c("0","10", "20", "30", "40", "50", "60", "70")) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  labs(fill = "Species") +
  ylab(label = "Number of scaffolds") +
  scale_fill_manual(values=c(wes_palette("Cavalcanti1")[2],wes_palette("Cavalcanti1")[3]))


# Obtain cumulative values of scaffold length
lei_i <-sort(lei$V2, decreasing = F)
fas_i <- sort(fas$length_fas, decreasing = F)
lei_i <- cbind(c(1:10433), lei_i)
fas_i <- cbind(c(1:2321), fas_i)
c <- c()
x = 0
for (i in 1:nrow(lei_i)){
  x = x + lei_i[i,2]
  c <- c(c,x)
}

lei_cumul <- as.data.frame(c)
colnames(lei_cumul) <- c("T.leiboldiana")
lei_cumul <- stack(lei_cumul)
lei_cumul$values <- sort(lei_cumul$values, decreasing = T)

c <- c()
x = 0
for (i in 1:nrow(fas_i)){
  x = x + fas_i[i,2]
  c <- c(c,x)
}
fas_cumul <- as.data.frame(c)
colnames(fas_cumul) <- c("T.fasciculata")
fas_cumul <- stack(fas_cumul)
fas_cumul$values <- sort(fas_cumul$values, decreasing = T)
cumulative_both <- rbind(fas_cumul, lei_cumul)
cumulative_both <- cbind(c(1:2321,1:10433), cumulative_both)
colnames(cumulative_both) <- c("scaffold", "length", "species")

ggplot(cumulative_both, aes(x = scaffold, y=length, group=species)) + geom_line()
cumulative_50 <- cumulative_both[c(1:50,2322:2371),]
ggplot(cumulative_50, aes(x = scaffold, y=length, group=species)) + geom_line()


# Plot cumulative distribution
cumul <- read.table("cumlength.csv", sep = ";", header = T)
ggplot(cumul, aes(x = scaffold, y=cumulative, group=species, colour=species)) + 
  geom_line() +
  scale_y_continuous(name = "Cumulative scaffold size (Mb)",
                     breaks=seq(0,1000000000,100000000),
                     labels = c("0","100", "200", "300", "400", "500", "600", "700", "800", "900", "1000")) +
  xlab(label="Number of scaffolds") +
  scale_color_manual(values=c(wes_palette("Cavalcanti1")[3],wes_palette("Cavalcanti1")[2]))

                   