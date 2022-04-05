setwd("bio-info_phd/DGE_Tfas-Tlei/")

d <- read.csv("go_enrichment.csv", header = T, sep = ";")

#install.packages("ggplot2")
#install.packages("wesanderson")
library(ggplot2)
library(wesanderson)
pal1 <- wes_palette("Darjeeling1")
pal2 <- wes_palette("Darjeeling2")
pal1[3]

ggplot(data = d, aes(x = species, y = carbon_fixation, fill = time)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.3) + 
  labs(title = "Carbon fixation", x = "Species", y = "Number of upregulated genes") +
  theme(plot.title = element_text(hjust = 0.5, family = "serif", face = "italic"),
        axis.title = element_text(family = "serif", vjust = 0.4),
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif")) +
  scale_fill_manual(values=c(pal1[3], pal2[2])) +
  coord_flip()

ggplot(data = d, aes(x = species, y = water_deprivation, fill = time)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.3) + 
  labs(title = "Response to water deprivation", x = "Species", y = "Number of upregulated genes") +
  theme(plot.title = element_text(hjust = 0.5, family = "serif", face = "italic"),
        axis.title = element_text(family = "serif", vjust = 0.4),
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif")) +
  scale_fill_manual(values=c(pal1[3], pal2[2])) +
  coord_flip()

ggplot(data = d, aes(x = species, y = photosynthesis_regulation, fill = time)) + 
  geom_bar(stat = "identity", position = position_dodge(), width = 0.3) + 
  labs(title = "Response to water deprivation", x = "Species", y = "Number of upregulated genes") +
  theme(plot.title = element_text(hjust = 0.5, family = "serif", face = "italic"),
        axis.title = element_text(family = "serif", vjust = 0.4),
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif")) +
  scale_fill_manual(values=c(pal1[3], pal2[2])) +
  coord_flip()

ggplot(data = d, aes(x = species, y = circadian_rhythm, fill = time)) + 
  geom_bar(stat = "identity", position = position_dodge())

