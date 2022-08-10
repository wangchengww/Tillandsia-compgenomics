#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2", "reshape2", "stringr","grid", "gridExtra", "matrixStats", "cowplot")

# Extract legend to plot legend for multiple plots
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Obtain arguments from command line. 1 = counts in T. fasciculata, 2 = counts in
# T. leiboldiana, 3 = gene list for orthogroup, 4 = color Tfas, 5 = color Tlei
args <- commandArgs(trailingOnly = TRUE)
print(args)
