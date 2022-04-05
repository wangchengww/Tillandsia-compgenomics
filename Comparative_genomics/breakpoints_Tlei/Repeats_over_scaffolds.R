setwd("/Users/clara/bio-info_phd/Comparative_genomics/breakpoints_Tlei/")

#install.packages("seqinr")
library("seqinr")
library("wesanderson")
Sc612 <- read.table(file = "")
Sc612a_seq <- Sc612[[3]]
Sc612b_seq <- Sc612[[2]]
Sc612c_seq <- Sc612[[1]]
Sc612d_seq <- Sc612[[4]]
Sc612e_seq <- Sc612[[5]]
Sc612f_seq <- Sc612[[6]]

# Scaffold_612a
out_a=NULL
starts_win <- seq(1, length(Sc612a_seq)-150000, by = 150000)
n <- length(starts_win)
for (i in 1:n) {
  subsetwindow <- Sc612a_seq[starts_win[i]:(starts_win[i]+149999)]
  a<- sum(lengths(regmatches(subsetwindow, gregexpr("a", subsetwindow))))
  c<- sum(lengths(regmatches(subsetwindow, gregexpr("c", subsetwindow))))
  g<- sum(lengths(regmatches(subsetwindow, gregexpr("g", subsetwindow))))
  t<- sum(lengths(regmatches(subsetwindow, gregexpr("t", subsetwindow))))
  masked_content <- sum(a,c,g,t)/150000
  out_a=rbind(out_a,cbind(starts_win[i],masked_content))
}

out_a=as.data.frame(out_a)
colnames(out_a)=c("position","masked_content")

plot(100*out_a$masked_content~out_a$position,type="l",lwd=4,lty=1,ylim=c(0,100),xlab="position",
     ylab="Proportion of masked bases", col=wes_palette("Darjeeling2")[2])+
  title(main = "Density of repetitive elements across 100 kb windows - Scaffold 612a")

#Scaffold_612b
out_b=NULL
starts_win <- seq(1, length(Sc612b_seq)-150000, by = 150000)
n <- length(starts_win)
for (i in 1:n) {
  subsetwindow <- Sc612b_seq[starts_win[i]:(starts_win[i]+149999)]
  a<- sum(lengths(regmatches(subsetwindow, gregexpr("a", subsetwindow))))
  c<- sum(lengths(regmatches(subsetwindow, gregexpr("c", subsetwindow))))
  g<- sum(lengths(regmatches(subsetwindow, gregexpr("g", subsetwindow))))
  t<- sum(lengths(regmatches(subsetwindow, gregexpr("t", subsetwindow))))
  masked_content <- sum(a,c,g,t)/150000
  out_b=rbind(out_b,cbind(starts_win[i],masked_content))
}

out_b=as.data.frame(out_b)
colnames(out_b)=c("position","masked_content")

plot(100*out_b$masked_content~out_b$position,type="l",lwd=4,lty=1,ylim=c(0,100),xlab="position",
     ylab="Proportion of masked bases", col=wes_palette("Darjeeling2")[2])+
  title(main = "Density of repetitive elements across 100 kb windows - Scaffold 612b")

# Scaffold_612c
out_c=NULL
starts_win <- seq(1, length(Sc612c_seq)-150000, by = 150000)
n <- length(starts_win)
for (i in 1:n) {
  subsetwindow <- Sc612c_seq[starts_win[i]:(starts_win[i]+149999)]
  a<- sum(lengths(regmatches(subsetwindow, gregexpr("a", subsetwindow))))
  c<- sum(lengths(regmatches(subsetwindow, gregexpr("c", subsetwindow))))
  g<- sum(lengths(regmatches(subsetwindow, gregexpr("g", subsetwindow))))
  t<- sum(lengths(regmatches(subsetwindow, gregexpr("t", subsetwindow))))
  masked_content <- sum(a,c,g,t)/150000
  out_c=rbind(out_c,cbind(starts_win[i],masked_content))
}

out_c=as.data.frame(out_c)
colnames(out_c)=c("position","masked_content")

plot(100*out_c$masked_content~out_c$position,type="l",lwd=4,lty=1,ylim=c(0,100),xlab="position",
     ylab="Proportion of masked bases", col=wes_palette("Darjeeling2")[2])+
  title(main = "Density of repetitive elements across 100 kb windows - Scaffold 612c")

# Scaffold_612d
out_d=NULL
starts_win <- seq(1, length(Sc612d_seq)-150000, by = 150000)
n <- length(starts_win)
for (i in 1:n) {
  subsetwindow <- Sc612d_seq[starts_win[i]:(starts_win[i]+149999)]
  a<- sum(lengths(regmatches(subsetwindow, gregexpr("a", subsetwindow))))
  c<- sum(lengths(regmatches(subsetwindow, gregexpr("c", subsetwindow))))
  g<- sum(lengths(regmatches(subsetwindow, gregexpr("g", subsetwindow))))
  t<- sum(lengths(regmatches(subsetwindow, gregexpr("t", subsetwindow))))
  masked_content <- sum(a,c,g,t)/150000
  out_d=rbind(out_d,cbind(starts_win[i],masked_content))
}

out_d=as.data.frame(out_d)
colnames(out_d)=c("position","masked_content")

plot(100*out_d$masked_content~out_d$position,type="l",lwd=4,lty=1,ylim=c(0,100),xlab="position",
     ylab="Proportion of masked bases", col=wes_palette("Darjeeling2")[2])+
  title(main = "Density of repetitive elements across 100 kb windows - Scaffold 612d")

# Scaffold_612e
out_e=NULL
starts_win <- seq(1, length(Sc612e_seq)-150000, by = 150000)
n <- length(starts_win)
for (i in 1:n) {
  subsetwindow <- Sc612e_seq[starts_win[i]:(starts_win[i]+149999)]
  a<- sum(lengths(regmatches(subsetwindow, gregexpr("a", subsetwindow))))
  c<- sum(lengths(regmatches(subsetwindow, gregexpr("c", subsetwindow))))
  g<- sum(lengths(regmatches(subsetwindow, gregexpr("g", subsetwindow))))
  t<- sum(lengths(regmatches(subsetwindow, gregexpr("t", subsetwindow))))
  masked_content <- sum(a,c,g,t)/150000
  out_e=rbind(out_e,cbind(starts_win[i],masked_content))
}

out_e=as.data.frame(out_e)
colnames(out_e)=c("position","masked_content")

plot(100*out_e$masked_content~out_e$position,type="l",lwd=4,lty=1,ylim=c(0,100),xlab="position",
     ylab="Proportion of masked bases", col=wes_palette("Darjeeling2")[2])+
  title(main = "Density of repetitive elements across 100 kb windows - Scaffold 612e")

# Scaffold_612f
out_f=NULL
starts_win <- seq(1, length(Sc612f_seq)-150000, by = 150000)
n <- length(starts_win)
for (i in 1:n) {
  subsetwindow <- Sc612f_seq[starts_win[i]:(starts_win[i]+149999)]
  a<- sum(lengths(regmatches(subsetwindow, gregexpr("a", subsetwindow))))
  c<- sum(lengths(regmatches(subsetwindow, gregexpr("c", subsetwindow))))
  g<- sum(lengths(regmatches(subsetwindow, gregexpr("g", subsetwindow))))
  t<- sum(lengths(regmatches(subsetwindow, gregexpr("t", subsetwindow))))
  masked_content <- sum(a,c,g,t)/150000
  out_f=rbind(out_f,cbind(starts_win[i],masked_content))
}

out_f=as.data.frame(out_f)
colnames(out_f)=c("position","masked_content")

plot(100*out_f$masked_content~out_f$position,type="l",lwd=4,lty=1,ylim=c(0,100),xlab="position",
     ylab="Proportion of masked bases", col=wes_palette("Darjeeling2")[2])+
  title(main = "Density of repetitive elements across 100 kb windows - Scaffold 612f")

