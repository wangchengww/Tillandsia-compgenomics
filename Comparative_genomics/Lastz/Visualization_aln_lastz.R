### Visualizing LastZ alignments
setwd("/Users/clara/bio-info_phd/Comparative_genomics/Lastz/")

### ALL CHROMOSOMES ###
anchors=read.table("Tlei_vs_Tfas_allchrom_lastz.filtered.0.9uniqfilter.coord",header=T)
require(reshape)

is.odd <- function(x) x %% 2 != 0

pdf("Tlei_vs_Tfas_allchrom_0.9uniqfilter2.pdf",10,10)
chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25")
library(RColorBrewer)
m <- 26
chromcolors <- colorRampPalette(brewer.pal(8, "Set1"))(m)
for (k in chr){
  if (k == "1") {
    limitchr=32162447
  } 
  else if (k == "2") {
    limitchr=31698969
  }
  else if (k == "3") {
    limitchr=31153634
  }
  else if (k == "4") {
    limitchr=30972320
  }
  else if (k == "5") {
    limitchr=30145683
  }
  else if (k == "6") {
    limitchr=29567185
  }
  else if (k == "7") {
    limitchr=27230796
  }
  else if (k == "8") {
    limitchr=27210704
  }
  else if (k == "9") {
    limitchr=26509606
  }
  else if (k == "10") {
    limitchr=26021936
  }
  else if (k == "11") {
    limitchr=25842354
  }
  else if (k == "12") {
    limitchr=25341332
  }
  else if (k == "13") {
    limitchr=24881768
  }
  else if (k == "14") {
    limitchr=24649150
  }
  else if (k == "15") {
    limitchr=23784123
  }
  else if (k == "16") {
    limitchr=23642332
  }
  else if (k == "17") {
    limitchr=21704838
  }
  else if (k == "18") {
    limitchr=21234157
  }
  else if (k == "19") {
    limitchr=19440187
  }
  else if (k == "20") {
    limitchr=19367790
  }
  else if (k == "21") {
    limitchr=18018563
  }
  else if (k == "22") {
    limitchr=17107477
  }
  else if (k == "23") {
    limitchr=16741216
  }
  else if (k == "24") {
    limitchr=13520883
  }
  else if (k == "25") {
    limitchr=12332728
  }
  chromo=paste("Tfas_chr",k,sep="")
  chromoshort=gsub("Tfas_","",gsub(",", "", chromo))
  print(chromo)
  plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(0, 1), ylim=c(0, 1),main=paste("\nTfas ",chromoshort," (",limitchr," bp)",sep=""))

  # add a scale
  xmin=0 # minor marks
  xbig=0 #main marks
  if (limitchr > 1000000){
    while (xmin < limitchr) {
      segments(x0=xmin/limitchr,y0=1.02,x1=xmin/limitchr,y1=1.015,col="darkgrey")
      xmin = xmin+1000000
    }
    while (xbig < limitchr) {
      segments(x0=xbig/limitchr,y0=1.02,x1=xbig/limitchr,y1=1.01,lwd=2)
      xbig = xbig+10000000
    }
  }
  segments(x0=0,y0=1.02,x1=limitchr/limitchr,y1=1.02,col="black",lwd=3)
  segments(x0=0,y0=0.5,x1=limitchr/limitchr,y1=0.5,col="black",lwd=1,lty=2)
  ### for each hits
  j=0 
  
  for (i in 1:nrow(anchors)){
    if (anchors[i,1] == chromo){
      j=j+1 # found hits
      startchr=anchors[i,2]
      endchr=anchors[i,2]+anchors[i,3]
      scaff=anchors[i,5]
      print(scaff)
      #scaff=gsub( "scaffold", "", scaff ) # replace scaffold
      scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
      scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
      scaff_n=as.numeric(gsub("chr", "", scaff))
      
      #print(scaff)
      #print(endchr)
      if (is.odd(j)=="TRUE"){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,lwd=3,col=chromcolors[scaff_n])
        segments(x0=(((endchr/limitchr)+(startchr/limitchr))/2),y0=0.5+decalline,x1=(((endchr/limitchr)+(startchr/limitchr))/2),y1=0.55+decal,col=chromcolors[scaff_n],lty=3)
        text((((endchr/limitchr)+(startchr/limitchr))/2),0.56+decal,labels=scaff[1],cex=0.7,col=chromcolors[scaff_n])
      }
      else {
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,col=chromcolors[scaff_n],lwd=3)
      segments(x0=(((endchr/limitchr)+(startchr/limitchr))/2),y0=0.5+decalline,x1=(((endchr/limitchr)+(startchr/limitchr))/2),y1=0.45-decal,col=chromcolors[scaff_n],lty=3)
      text((((endchr/limitchr)+(startchr/limitchr))/2),0.44-decal,labels=scaff[1],cex=0.7,col=chromcolors[scaff_n])
      }
      #print(anchors)
      #print(anchors[i,])
    }
  }
}
dev.off()

##################################
## Make PDF just for Tlei_chr14 ##
##################################
tlei14 <- read.table("Tfas_vs_Tlei_chr14_LASTZ_alignment0.9uniqfilter.coord", header = T)
pdf("Tlei_chr14_0.9uniqfilter_greens.pdf",10,10)
k = 14
limitchr=38474115
chromo=paste("Tlei_chr",k,sep="")
chromoshort=gsub("Tlei_","",gsub(",", "", chromo))
print(chromo)
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(0, 1.1), ylim=c(0, .7),main=paste("\nTlei ",chromoshort," (",limitchr," bp)",sep=""))

library(RColorBrewer)
m <- 24
chromcolors <- colorRampPalette(brewer.pal(8, "Greys"))(m)

# add a scale
xmin=0 # minor marks
xbig=0 #main marks
if (limitchr > 1000000){
  while (xmin < limitchr) {
    segments(x0=xmin/limitchr,y0=0.62,x1=xmin/limitchr,y1=.615,col="darkgrey")
    xmin = xmin+1000000
  }
  while (xbig < limitchr) {
    segments(x0=xbig/limitchr,y0=.62,x1=xbig/limitchr,y1=.61,lwd=2)
    text(x = (xbig/limitchr), y = 0.6, labels = xbig, col="black")
    xbig = xbig+10000000
  }
}
segments(x0=0,y0=.62,x1=limitchr/limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr/limitchr,y1=0.5,col="black",lwd=1,lty=2)

### for each hits
j=0 

for (i in 1:nrow(tlei14)){
  if (tlei14[i,1] == chromo){
    j=j+1 # found hits
    startchr=tlei14[i,2]
    endchr=tlei14[i,2]+tlei14[i,3]
    scaff=tlei14[i,5]
    print(scaff)
    scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
    scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
    scaff_n=as.numeric(gsub("chr", "", scaff))
    
    if (is.odd(j)=="TRUE"){
      
      if (scaff_n == 17){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,lwd=3,col="maroon3")

      } else if (scaff_n == 25){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,lwd=3,col="orange2")

      } else {
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,lwd=3,col=chromcolors[scaff_n])
      }
    }
    else {
      if (scaff_n == 17){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,col="maroon3",lwd=3)
        
      } else if (scaff_n == 25){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,col="orange2",lwd=3)
      } else {
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,col=chromcolors[scaff_n],lwd=3)
      }
    }
    #print(anchors)
    #print(anchors[i,])
  }
}
text(x = 1.05, y = 0.55, label = "Tfas_chr17", col="black", cex = .6)
text(x = 1.05, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 1.05, y = 0.45, label = "Tfas_chr25", col="black", cex = .6)
dev.off()

### Zoom in to the breakpoint ###
zoom <- tlei14[tlei14$Qpos > 28400000 & tlei14$Qpos < 28700000,]
pdf("Tlei_chr14_zoom_300000_greys.pdf",10,10)
# add a scale
limitchr = 28.7
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(28.4, 28.75), ylim=c(0, .7),main=paste("\nTlei chr14 (28,4 - 28,7 MB)",sep=""))
xmin=28.4 # minor marks
xbig=28.4 #main marks
while (xmin < limitchr) {
  segments(x0=xmin,y0=0.62,x1=xmin,y1=.615,col="darkgrey")
  xmin = xmin+0.01
}
while (xbig < limitchr) {
  segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
  text(x = (xbig), y = 0.6, labels = xbig, col="black")
  xbig = xbig+0.1
}
segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
text(x = (xbig), y = 0.6, labels = xbig, col="black")
segments(x0=0,y0=.62,x1=limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr,y1=0.5,col="black",lwd=1,lty=2)

j=0 
for (i in 1:nrow(zoom)){
  j=j+1 # found hits
  startchr=zoom[i,2]/1000000
  endchr=(zoom[i,2]+zoom[i,3])/1000000
  scaff=zoom[i,5]
  print(scaff)
  scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
  scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
  scaff_n=as.numeric(gsub("chr", "", scaff))
  
  if (is.odd(j)=="TRUE"){
    
    if (scaff_n == 17){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,lwd=3,col="maroon3")
      
    } else if (scaff_n == 25){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,lwd=3,col="orange2")
      
    } else {
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,lwd=3,col=chromcolors[scaff_n])
    }
  }
  else {
    if (scaff_n == 17){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,col="maroon3",lwd=3)
      
    } else if (scaff_n == 25){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,col="orange2",lwd=3)
    } else {
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,col=chromcolors[scaff_n],lwd=3)
    }
  }
}
text(x = 28.72, y = 0.55, label = "Tfas_chr17", col="black", cex = .6)
text(x = 28.72, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 28.72, y = 0.45, label = "Tfas_chr25", col="black", cex = .6)
dev.off()

############################################
## Make PDF just for chr2 & chr13 in Tfas ##
############################################

pdf("Tfas_chr2-13_0.9uniqfilter2.pdf",10,10)
k = 2
limitchr=31698969
chromo=paste("Tfas_chr",k,sep="")
chromoshort=gsub("Tfas_","",gsub(",", "", chromo))
print(chromo)
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(0, 1.1), ylim=c(0, .7),main=paste("\nTfas ",chromoshort," (",limitchr," bp)",sep=""))

# add a scale
xmin=0 # minor marks
xbig=0 #main marks
if (limitchr > 1000000){
  while (xmin < limitchr) {
    segments(x0=xmin/limitchr,y0=0.62,x1=xmin/limitchr,y1=.615,col="darkgrey")
    xmin = xmin+1000000
  }
  while (xbig < limitchr) {
    segments(x0=xbig/limitchr,y0=.62,x1=xbig/limitchr,y1=.61,lwd=2)
    text(x = (xbig/limitchr), y = 0.6, labels = xbig, col="black")
    xbig = xbig+10000000
  }
}
segments(x0=0,y0=.62,x1=limitchr/limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr/limitchr,y1=0.5,col="black",lwd=1,lty=2)

### for each hits
j=0 
for (i in 1:nrow(anchors)){
  if (anchors[i,1] == chromo){
    j=j+1 # found hits
    startchr=anchors[i,2]
    endchr=anchors[i,2]+anchors[i,3]
    scaff=anchors[i,5]
    print(scaff)
    scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
    scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
    scaff_n=as.numeric(gsub("chr", "", scaff))
    
    if (is.odd(j)=="TRUE"){
      
      if (scaff_n == 19){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,lwd=3,col="turquoise3")
        
      } else if (scaff_n == 3){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,lwd=3,col="orange2")
        
      } else {
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,lwd=3,col="grey50")
      }
    }
    else {
      if (scaff_n == 19){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,col="turquoise3",lwd=3)
        
      } else if (scaff_n == 3){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,col="orange2",lwd=3)
      } else {
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,col="grey50",lwd=3)
      }
    }
    #print(anchors)
    #print(anchors[i,])
  }
}
text(x = 1.05, y = 0.55, label = "Tlei_chr19", col="black", cex = .6)
text(x = 1.05, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 1.05, y = 0.45, label = "Tlei_chr3", col="black", cex = .6)

k = 13
limitchr=24881768
chromo=paste("Tfas_chr",k,sep="")
chromoshort=gsub("Tfas_","",gsub(",", "", chromo))
print(chromo)
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(0, 1.1), ylim=c(0, .7),main=paste("\nTfas ",chromoshort," (",limitchr," bp)",sep=""))
tfas13 <- anchors[anchors$Qname=="Tfas_chr13",]
# add a scale
xmin=0 # minor marks
xbig=0 #main marks
if (limitchr > 1000000){
  while (xmin < limitchr) {
    segments(x0=xmin/limitchr,y0=0.62,x1=xmin/limitchr,y1=.615,col="darkgrey")
    xmin = xmin+1000000
  }
  while (xbig < limitchr) {
    segments(x0=xbig/limitchr,y0=.62,x1=xbig/limitchr,y1=.61,lwd=2)
    text(x = (xbig/limitchr), y = 0.6, labels = xbig, col="black")
    xbig = xbig+10000000
  }
}
segments(x0=0,y0=.62,x1=limitchr/limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr/limitchr,y1=0.5,col="black",lwd=1,lty=2)

### for each hits
j=0 
for (i in 1:nrow(anchors)){
  if (anchors[i,1] == chromo){
    j=j+1 # found hits
    startchr=anchors[i,2]
    endchr=anchors[i,2]+anchors[i,3]
    scaff=anchors[i,5]
    print(scaff)
    scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
    scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
    scaff_n=as.numeric(gsub("chr", "", scaff))
    
    if (is.odd(j)=="TRUE"){
      
      if (scaff_n == 19){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,lwd=3,col="turquoise3")
        
      } else if (scaff_n == 3){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,lwd=3,col="orange2")
        
      } else {
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,lwd=3,col="grey50")
      }
    }
    else {
      if (scaff_n == 19){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,col="turquoise3",lwd=3)
        
      } else if (scaff_n == 3){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,col="orange2",lwd=3)
      } else {
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,col="grey50",lwd=3)
      }
    }
    #print(anchors)
    #print(anchors[i,])
  }
}
text(x = 1.05, y = 0.55, label = "Tlei_chr19", col="black", cex = .6)
text(x = 1.05, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 1.05, y = 0.45, label = "Tlei_chr3", col="black", cex = .6)
dev.off()

###
# Zoom in on centromere chromosome 2 
zoomchr2 <- anchors[anchors$Qname == "Tfas_chr2" & anchors$Qpos > 15000000 & anchors$Qpos < 23000000,]
pdf("Tfas_chr2_zoom_8MB.pdf",10,10)
# add a scale
limitchr = 23
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(15, 23.5), ylim=c(0, .7),main=paste("\nTfas chr2 (15 - 23 MB)",sep=""))
xmin=15 # minor marks
xbig=15 #main marks
while (xmin < limitchr) {
  segments(x0=xmin,y0=0.62,x1=xmin,y1=.615,col="darkgrey")
  xmin = xmin+1
}
while (xbig < limitchr) {
  segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
  text(x = (xbig), y = 0.6, labels = xbig, col="black")
  xbig = xbig+4
}
segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
text(x = (xbig), y = 0.6, labels = xbig, col="black")
segments(x0=0,y0=.62,x1=limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr,y1=0.5,col="black",lwd=1,lty=2)

j=0 
for (i in 1:nrow(zoomchr2)){
  j=j+1 # found hits
  startchr=zoomchr2[i,2]/1000000
  endchr=(zoomchr2[i,2]+zoomchr2[i,3])/1000000
  scaff=zoomchr2[i,5]
  print(scaff)
  scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
  scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
  scaff_n=as.numeric(gsub("chr", "", scaff))
  
  if (is.odd(j)=="TRUE"){
    
    if (scaff_n == 19){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,lwd=3,col="turquoise3")
      
    } else if (scaff_n == 3){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,lwd=3,col="orange2")
      
    } else {
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,lwd=3,col="grey50")
    }
  }
  else {
    if (scaff_n == 19){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,col="turquoise3",lwd=3)
      
    } else if (scaff_n == 3){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,col="orange2",lwd=3)
    } else {
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,col="grey50",lwd=3)
    }
  }
}
text(x = 23.5, y = 0.55, label = "Tlei_chr19", col="black", cex = .6)
text(x = 23.5, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 23.5, y = 0.45, label = "Tlei_chr3", col="black", cex = .6)
dev.off()

###
# Zoom in on breakpoints chr13
###
# Breakpoint 1 
zoomchr13 <- anchors[anchors$Qname == "Tfas_chr13" & anchors$Qpos > 1000000 & anchors$Qpos < 1500000,]
pdf("Tfas_chr13_zoom_500MB_breakpoint1.pdf",10,10)
# add a scale
limitchr = 1.5
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(1, 1.55), ylim=c(0, .7),main=paste("\nTfas chr13 (1 - 1.5 MB)",sep=""))
xmin=1 # minor marks
xbig=1 #main marks
while (xmin < limitchr) {
  segments(x0=xmin,y0=0.62,x1=xmin,y1=.615,col="darkgrey")
  xmin = xmin+0.01
}
while (xbig < limitchr) {
  segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
  text(x = (xbig), y = 0.6, labels = xbig, col="black")
  xbig = xbig+0.1
}
segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
text(x = (xbig), y = 0.6, labels = xbig, col="black")
segments(x0=0,y0=.62,x1=limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr,y1=0.5,col="black",lwd=1,lty=2)

j=0 
for (i in 1:nrow(zoomchr13)){
  j=j+1 # found hits
  startchr=zoomchr13[i,2]/1000000
  endchr=(zoomchr13[i,2]+zoomchr13[i,3])/1000000
  scaff=zoomchr13[i,5]
  print(scaff)
  scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
  scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
  scaff_n=as.numeric(gsub("chr", "", scaff))
  
  if (is.odd(j)=="TRUE"){
    
    if (scaff_n == 19){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,lwd=3,col="turquoise3")
      
    } else if (scaff_n == 3){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,lwd=3,col="orange2")
      
    } else {
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,lwd=3,col="grey50")
    }
  }
  else {
    if (scaff_n == 19){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,col="turquoise3",lwd=3)
      
    } else if (scaff_n == 3){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,col="orange2",lwd=3)
    } else {
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,col="grey50",lwd=3)
    }
  }
}
text(x = 1.55, y = 0.55, label = "Tlei_chr19", col="black", cex = .6)
text(x = 1.55, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 1.55, y = 0.45, label = "Tlei_chr3", col="black", cex = .6)
dev.off()

###
# Breakpoint 2
zoomchr13b <- anchors[anchors$Qname == "Tfas_chr13" & anchors$Qpos > 4100000 & anchors$Qpos < 4200000,]
pdf("Tfas_chr13_zoom_100KB_breakpoint2.pdf",10,10)
# add a scale
limitchr = 4.2
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(4.1, 4.21), ylim=c(0, .7),main=paste("\nTfas chr13 (4.1 - 4.2 MB)",sep=""))
xmin=4.1 # minor marks
xbig=4.1 #main marks
while (xmin < limitchr) {
  segments(x0=xmin,y0=0.62,x1=xmin,y1=.615,col="darkgrey")
  xmin = xmin+0.01
}
while (xbig < limitchr) {
  segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
  text(x = (xbig), y = 0.6, labels = xbig, col="black")
  xbig = xbig+0.1
}
segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
text(x = (xbig), y = 0.6, labels = xbig, col="black")
segments(x0=0,y0=.62,x1=limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr,y1=0.5,col="black",lwd=1,lty=2)

j=0 
for (i in 1:nrow(zoomchr13b)){
  j=j+1 # found hits
  startchr=zoomchr13b[i,2]/1000000
  endchr=(zoomchr13b[i,2]+zoomchr13b[i,3])/1000000
  scaff=zoomchr13b[i,5]
  print(scaff)
  scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
  scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
  scaff_n=as.numeric(gsub("chr", "", scaff))
  
  if (is.odd(j)=="TRUE"){
    
    if (scaff_n == 19){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,lwd=3,col="turquoise3")
      
    } else if (scaff_n == 3){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,lwd=3,col="orange2")
      
    } else {
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,lwd=3,col="grey50")
    }
  }
  else {
    if (scaff_n == 19){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,col="turquoise3",lwd=3)
      
    } else if (scaff_n == 3){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,col="orange2",lwd=3)
    } else {
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,col="grey50",lwd=3)
    }
  }
}
text(x = 4.21, y = 0.55, label = "Tlei_chr19", col="black", cex = .6)
text(x = 4.21, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 4.21, y = 0.45, label = "Tlei_chr3", col="black", cex = .6)
dev.off()

#############################################
## Make PDF just for chr10 & chr24 in Tfas ##
#############################################

pdf("Tfas_chr10-24_0.9uniqfilter.pdf",10,10)
k = 10
limitchr=26021936
chromo=paste("Tfas_chr",k,sep="")
chromoshort=gsub("Tfas_","",gsub(",", "", chromo))
print(chromo)
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(0, 1.1), ylim=c(0, .7),main=paste("\nTfas ",chromoshort," (",limitchr," bp)",sep=""))

# add a scale
xmin=0 # minor marks
xbig=0 #main marks
if (limitchr > 1000000){
  while (xmin < limitchr) {
    segments(x0=xmin/limitchr,y0=0.62,x1=xmin/limitchr,y1=.615,col="darkgrey")
    xmin = xmin+1000000
  }
  while (xbig < limitchr) {
    segments(x0=xbig/limitchr,y0=.62,x1=xbig/limitchr,y1=.61,lwd=2)
    text(x = (xbig/limitchr), y = 0.6, labels = xbig, col="black")
    xbig = xbig+10000000
  }
}
segments(x0=0,y0=.62,x1=limitchr/limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr/limitchr,y1=0.5,col="black",lwd=1,lty=2)

### for each hits
j=0 
for (i in 1:nrow(anchors)){
  if (anchors[i,1] == chromo){
    j=j+1 # found hits
    startchr=anchors[i,2]
    endchr=anchors[i,2]+anchors[i,3]
    scaff=anchors[i,5]
    print(scaff)
    scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
    scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
    scaff_n=as.numeric(gsub("chr", "", scaff))
    
    if (is.odd(j)=="TRUE"){
      
      if (scaff_n == 13){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,lwd=3,col="springgreen3")
        
      } else if (scaff_n == 23){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,lwd=3,col="tomato2")
        
      } else {
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,lwd=3,col="grey50")
      }
    }
    else {
      if (scaff_n == 13){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,col="springgreen3",lwd=3)
        
      } else if (scaff_n == 23){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,col="tomato2",lwd=3)
      } else {
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,col="grey50",lwd=3)
      }
    }
    #print(anchors)
    #print(anchors[i,])
  }
}
text(x = 1.05, y = 0.55, label = "Tlei_chr19", col="black", cex = .6)
text(x = 1.05, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 1.05, y = 0.45, label = "Tlei_chr3", col="black", cex = .6)

k = 24
limitchr=13520883
chromo=paste("Tfas_chr",k,sep="")
chromoshort=gsub("Tfas_","",gsub(",", "", chromo))
print(chromo)
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(0, 1.1), ylim=c(0, .7),main=paste("\nTfas ",chromoshort," (",limitchr," bp)",sep=""))

# add a scale
xmin=0 # minor marks
xbig=0 #main marks
if (limitchr > 1000000){
  while (xmin < limitchr) {
    segments(x0=xmin/limitchr,y0=0.62,x1=xmin/limitchr,y1=.615,col="darkgrey")
    xmin = xmin+1000000
  }
  while (xbig < limitchr) {
    segments(x0=xbig/limitchr,y0=.62,x1=xbig/limitchr,y1=.61,lwd=2)
    text(x = (xbig/limitchr), y = 0.6, labels = xbig, col="black")
    xbig = xbig+10000000
  }
}
segments(x0=0,y0=.62,x1=limitchr/limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr/limitchr,y1=0.5,col="black",lwd=1,lty=2)

### for each hits
j=0 
for (i in 1:nrow(anchors)){
  if (anchors[i,1] == chromo){
    j=j+1 # found hits
    startchr=anchors[i,2]
    endchr=anchors[i,2]+anchors[i,3]
    scaff=anchors[i,5]
    print(scaff)
    scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
    scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
    scaff_n=as.numeric(gsub("chr", "", scaff))
    
    if (is.odd(j)=="TRUE"){
      
      if (scaff_n == 13){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,lwd=3,col="springgreen3")
        
      } else if (scaff_n == 23){
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,lwd=3,col="tomato2")
        
      } else {
        decalline=runif(1,0,0.015)
        decal=runif(1, decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,lwd=3,col="grey50")
      }
    }
    else {
      if (scaff_n == 13){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.55+decalline,x1=(endchr/limitchr),y1=0.55+decalline,col="springgreen3",lwd=3)
        
      } else if (scaff_n == 23){
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.45+decalline,x1=(endchr/limitchr),y1=0.45+decalline,col="tomato2",lwd=3)
      } else {
        decalline=runif(1,-0.015,0)
        decal=runif(1, -decalline, 0.4)
        segments(x0=(startchr/limitchr),y0=0.5+decalline,x1=(endchr/limitchr),y1=0.5+decalline,col="grey50",lwd=3)
      }
    }
    #print(anchors)
    #print(anchors[i,])
  }
}
text(x = 1.05, y = 0.55, label = "Tlei_chr13", col="black", cex = .6)
text(x = 1.05, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 1.05, y = 0.45, label = "Tlei_chr23", col="black", cex = .6)
dev.off()

###
# Zoom in to breakpoint chr 10
zoomchr10 <- anchors[anchors$Qname == "Tfas_chr10" & anchors$Qpos > 22800000 & anchors$Qpos < 23200000,]
pdf("Tfas_chr10_zoom_400KB.pdf",10,10)
# add a scale
limitchr = 23.2
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(22.8, 23.22), ylim=c(0, .7),main=paste("\nTfas chr10 (22.8 - 23.2 MB)",sep=""))
xmin=22.8 # minor marks
xbig=22.8 #main marks
while (xmin < limitchr) {
  segments(x0=xmin,y0=0.62,x1=xmin,y1=.615,col="darkgrey")
  xmin = xmin+0.01
}
while (xbig < limitchr) {
  segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
  text(x = (xbig), y = 0.6, labels = xbig, col="black")
  xbig = xbig+0.1
}
segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
text(x = (xbig), y = 0.6, labels = xbig, col="black")
segments(x0=0,y0=.62,x1=limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr,y1=0.5,col="black",lwd=1,lty=2)

j=0 
for (i in 1:nrow(zoomchr10)){
  j=j+1 # found hits
  startchr=zoomchr10[i,2]/1000000
  endchr=(zoomchr10[i,2]+zoomchr10[i,3])/1000000
  scaff=zoomchr10[i,5]
  print(scaff)
  scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
  scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
  scaff_n=as.numeric(gsub("chr", "", scaff))
  
  if (is.odd(j)=="TRUE"){
    
    if (scaff_n == 13){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,lwd=3,col="springgreen3")
      
    } else if (scaff_n == 23){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,lwd=3,col="tomato2")
      
    } else {
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,lwd=3,col="grey50")
    }
  }
  else {
    if (scaff_n == 13){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,col="springgreen3",lwd=3)
      
    } else if (scaff_n == 23){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,col="tomato2",lwd=3)
    } else {
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,col="grey50",lwd=3)
    }
  }
}
text(x = 23.22, y = 0.55, label = "Tlei_chr13", col="black", cex = .6)
text(x = 23.22, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 23.22, y = 0.45, label = "Tlei_chr23", col="black", cex = .6)
dev.off()

###
# Zoom in on chr24
zoomchr24 <- anchors[anchors$Qname == "Tfas_chr24" & anchors$Qpos > 4000000 & anchors$Qpos < 5000000,]
pdf("Tfas_chr24_zoom_1MB.pdf",10,10)
# add a scale
limitchr = 5
plot(NULL, type="n", axes=FALSE, xlab="",ylab="",xlim=c(4, 5.04), ylim=c(0, .7),main=paste("\nTfas chr24 (4 - 5 MB)",sep=""))
xmin=4 # minor marks
xbig=4 #main marks
while (xmin < limitchr) {
  segments(x0=xmin,y0=0.62,x1=xmin,y1=.615,col="darkgrey")
  xmin = xmin+0.01
}
while (xbig < limitchr) {
  segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
  text(x = (xbig), y = 0.6, labels = xbig, col="black")
  xbig = xbig+0.1
}
segments(x0=xbig,y0=.62,x1=xbig,y1=.61,lwd=2)
text(x = (xbig), y = 0.6, labels = xbig, col="black")
segments(x0=0,y0=.62,x1=limitchr,y1=.62,col="black",lwd=3)
segments(x0=0,y0=0.5,x1=limitchr,y1=0.5,col="black",lwd=1,lty=2)

j=0 
for (i in 1:nrow(zoomchr24)){
  j=j+1 # found hits
  startchr=zoomchr24[i,2]/1000000
  endchr=(zoomchr24[i,2]+zoomchr24[i,3])/1000000
  scaff=zoomchr24[i,5]
  print(scaff)
  scaff=strsplit(scaff,"_") # split the n°scaff & the size of the scaff (if exists)
  scaff=sapply(scaff,function(x) x[2]) # keep only the n°scaff
  scaff_n=as.numeric(gsub("chr", "", scaff))
  
  if (is.odd(j)=="TRUE"){
    
    if (scaff_n == 13){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,lwd=3,col="springgreen3")
      
    } else if (scaff_n == 23){
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,lwd=3,col="tomato2")
      
    } else {
      decalline=runif(1,0,0.015)
      decal=runif(1, decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,lwd=3,col="grey50")
    }
  }
  else {
    if (scaff_n == 13){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.55+decalline,x1=(endchr),y1=0.55+decalline,col="springgreen3",lwd=3)
      
    } else if (scaff_n == 23){
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.45+decalline,x1=(endchr),y1=0.45+decalline,col="tomato2",lwd=3)
    } else {
      decalline=runif(1,-0.015,0)
      decal=runif(1, -decalline, 0.4)
      segments(x0=(startchr),y0=0.5+decalline,x1=(endchr),y1=0.5+decalline,col="grey50",lwd=3)
    }
  }
}
text(x = 5.04, y = 0.55, label = "Tlei_chr13", col="black", cex = .6)
text(x = 5.04, y = 0.5, label = "Other chr", col="black", cex = .6)
text(x = 5.04, y = 0.45, label = "Tlei_chr23", col="black", cex = .6)
dev.off()
