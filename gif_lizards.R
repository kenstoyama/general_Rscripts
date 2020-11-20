setwd("~/Desktop/R/")
install.packages("ggplot2")
install.packages("ggfortify")
install.packages("animation")
install.packages("magick")
install.packages("gifski")
install.packages("ggpubr")
install.packages("ggforce")
library(ggpubr)
library(ggfortify)
library(gifski)
library(ggplot2)
library(animation)
library(magick)
library(gtools)
library(gifski)
library(ggforce)
data <- read.table("pruebaespecies.txt",header=T, na.strings = c("x", "NA"))
prueba.pca <- prcomp(data[1:154,3:16])
new <- prueba.pca$x[,1:2]
gen <- as.character(data$gen)
new <- as.data.frame(cbind(gen,new))
new$PC1 <- as.numeric(as.character(new$PC1))
new$PC2 <- as.numeric(as.character(new$PC2))

#### GIF METHOD 2 ####
for (i in 1:180){
  png(paste("test",i,".png"), units="in", width=8, height=4, res=200)
  degs <- seq(0,2*pi,by=(2*pi)/180)
  X <- -0.05+0.15*sin(degs[i]) # x coordinate
  Y <- -0.05+0.15*cos(degs[i]) # y coordinate
  a <- autoplot(prueba.pca,data=new,colour="gen",frame=T, xlim=c(-0.3,0.2)) + 
    geom_point(aes(x=X, y=Y), colour="firebrick",cex=4,pch=3, stroke=1.5) +
    coord_fixed() + guides(shape = guide_legend(override.aes = list(size = 0.5))) +
    theme(legend.title = element_blank(),plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          axis.title=element_text(size=10),legend.position="right")
  df <- data.frame()
  body.xa <- -0.2-((X*prueba.pca$rotation[4]+Y*prueba.pca$rotation[4,2]))
  body.xb <- -0.2-((X*prueba.pca$rotation[4]+Y*prueba.pca$rotation[4,2]))
  body.xc <- 0.2+((X*prueba.pca$rotation[4]+Y*prueba.pca$rotation[4,2]))
  body.xd <- 0.2+((X*prueba.pca$rotation[4]+Y*prueba.pca$rotation[4,2]))
  body.xe <- body.xd-0.05
  body.xf <- body.xa+0.05
  body=data.frame(x=c(body.xa,body.xb,body.xc,body.xd,body.xe,body.xf), y=c(-0.4,0.4,0.4,-0.4,-0.5,-0.5))
  head.xa <- -0.15-((X*prueba.pca$rotation[2]+Y*prueba.pca$rotation[1,2]))
  head.xc <- 0.15+((X*prueba.pca$rotation[2]+Y*prueba.pca$rotation[1,2]))
  head.yb <- 0.65+X*prueba.pca$rotation[1]+Y*prueba.pca$rotation[1,2]
  head=data.frame(x=c(head.xa,head.xa,0,head.xc,head.xc), y=c(0.4+0.025,0.5+0.025,head.yb,0.5+0.025,0.4+0.025))
  humerusL.xa <- body.xa-0.3-((X*prueba.pca$rotation[11]+Y*prueba.pca$rotation[11,2]))
  humerusL.xb <- body.xb-0.2-((X*prueba.pca$rotation[11]+Y*prueba.pca$rotation[11,2]))
  humerusL=data.frame(x=c(humerusL.xa,humerusL.xb,body.xa-0.025,body.xa-0.025), y=c(0.25,0.35,0.35,0.25))
  humerusR.xd <- body.xc+0.3+((X*prueba.pca$rotation[11]+Y*prueba.pca$rotation[11,2]))
  humerusR.xc <- body.xd+0.2+((X*prueba.pca$rotation[11]+Y*prueba.pca$rotation[11,2]))
  humerusR=data.frame(x=c(body.xd+0.025,body.xd+0.025,humerusR.xc,humerusR.xd), y=c(0.25,0.35,0.35,0.25))
  radiusL.xa <- humerusL.xa-0.01
  radiusL.xd <- humerusL.xb-0.01
  radiusL.yb <- 0.55+((X*prueba.pca$rotation[12]+Y*prueba.pca$rotation[12,2]))
  radiusL.yc <- 0.55+((X*prueba.pca$rotation[12]+Y*prueba.pca$rotation[12,2]))
  radiusL=data.frame(x=c(radiusL.xa,radiusL.xa,radiusL.xd,radiusL.xd), y=c(0.26,radiusL.yb,radiusL.yc,0.36))
  radiusR.xa <- humerusR.xc+0.01
  radiusR.xd <- humerusR.xd+0.01
  radiusR.yb <- 0.55+((X*prueba.pca$rotation[12]+Y*prueba.pca$rotation[12,2]))
  radiusR.yc <- 0.55+((X*prueba.pca$rotation[12]+Y*prueba.pca$rotation[12,2]))
  radiusR=data.frame(x=c(radiusR.xa,radiusR.xa,radiusR.xd,radiusR.xd), y=c(0.36,radiusR.yb,radiusR.yc,0.26))
  metacL.yb <- radiusL.yb+0.01+0.06+0.2*(X*prueba.pca$rotation[13]+Y*prueba.pca$rotation[13,2])
  metacL=data.frame(x=c(radiusL.xa-0.01,radiusL.xa-0.01,radiusL.xd+0.01,radiusL.xd+0.01), y=c(radiusL.yb+0.01,metacL.yb,metacL.yb,radiusL.yb+0.01))
  metacR.yb <- radiusR.yb+0.01+0.06+0.2*(X*prueba.pca$rotation[13]+Y*prueba.pca$rotation[13,2])
  metacR=data.frame(x=c(radiusR.xa-0.01,radiusR.xa-0.01,radiusR.xd+0.01,radiusR.xd+0.01), y=c(radiusR.yb+0.01,metacR.yb,metacR.yb,radiusR.yb+0.01))
  distlimbs <- (X*prueba.pca$rotation[6]+Y*prueba.pca$rotation[6,2])
  femurL.xa <- body.xa-0.3-((X*prueba.pca$rotation[7]+Y*prueba.pca$rotation[7,2]))
  femurL.xb <- body.xb-0.4-((X*prueba.pca$rotation[7]+Y*prueba.pca$rotation[7,2]))
  femurL=data.frame(x=c(femurL.xa,femurL.xb,body.xa-0.025,body.xa-0.025), y=c(-0.30-distlimbs,-0.20-distlimbs,-0.20-distlimbs,-0.30-distlimbs))
  femurR.xc <- body.xc+0.4+((X*prueba.pca$rotation[7]+Y*prueba.pca$rotation[7,2]))
  femurR.xd <- body.xd+0.3+((X*prueba.pca$rotation[7]+Y*prueba.pca$rotation[7,2]))
  femurR=data.frame(x=c(body.xd+0.025,body.xd+0.025,femurR.xc,femurR.xd), y=c(-0.30-distlimbs,-0.20-distlimbs,-0.20-distlimbs,-0.30-distlimbs))
  tibiaL.xa <- femurL.xb-0.01
  tibiaL.xb <- femurL.xb-0.01
  tibiaL.xc <- femurL.xa-0.01
  tibiaL.xd <- femurL.xa-0.01
  tibiaL.ya <- -0.21-0.3-(X*prueba.pca$rotation[8]+Y*prueba.pca$rotation[8,2])
  tibiaL=data.frame(x=c(tibiaL.xa,tibiaL.xb,tibiaL.xc,tibiaL.xd), y=c(tibiaL.ya-distlimbs,-0.21-distlimbs,-0.31-distlimbs,tibiaL.ya-distlimbs))
  tibiaR.xa <- femurR.xd+0.01
  tibiaR.xb <- femurR.xd+0.01
  tibiaR.xc <- femurR.xc+0.01
  tibiaR.xd <- femurR.xc+0.01
  tibiaR.ya <- -0.21-0.3-(X*prueba.pca$rotation[8]+Y*prueba.pca$rotation[8,2])
  tibiaR=data.frame(x=c(tibiaR.xa,tibiaR.xb,tibiaR.xc,tibiaR.xd), y=c(tibiaL.ya-distlimbs,-0.31-distlimbs,-0.21-distlimbs,tibiaL.ya-distlimbs))
  metatL.xa <- tibiaL.xa-0.25-(X*prueba.pca$rotation[9]+Y*prueba.pca$rotation[9,2])
  metatL=data.frame(x=c(metatL.xa,metatL.xa,tibiaL.xc,tibiaL.xc), y=c(tibiaL.ya-distlimbs-0.06,tibiaL.ya-distlimbs-0.01,tibiaL.ya-distlimbs-0.01,tibiaL.ya-distlimbs-0.06))
  metatR.xc <- tibiaR.xc+0.25+(X*prueba.pca$rotation[9]+Y*prueba.pca$rotation[9,2])
  metatR=data.frame(x=c(tibiaR.xa,tibiaR.xa,metatR.xc,metatR.xc), y=c(tibiaR.ya-distlimbs-0.06,tibiaL.ya-distlimbs-0.01,tibiaL.ya-distlimbs-0.01,tibiaL.ya-distlimbs-0.06))
  fingerHL1.x <- c(metatL.xa+0.02,metatL.xa-0.1-0.05*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHL1.y <- c(tibiaR.ya-distlimbs-0.05,tibiaR.ya-distlimbs-0.12-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHL1=data.frame(x=fingerHL1.x,y=fingerHL1.y)
  fingerHL2.x <- c(metatL.xa,metatL.xa-0.15-0.05*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHL2.y <- c(tibiaR.ya-distlimbs-0.025,tibiaR.ya-distlimbs-0.05-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHL2=data.frame(x=fingerHL2.x,y=fingerHL2.y)
  fingerHL3.x <- c(metatL.xa+0.1,metatL.xa-0.05*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHL3.y <- c(tibiaR.ya-distlimbs-0.025,tibiaR.ya-distlimbs-0.15-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHL3=data.frame(x=fingerHL3.x,y=fingerHL3.y)
  fingerHR1.x <- c(metatR.xc-0.02,metatR.xc+0.1-0.05*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHR1.y <- c(tibiaR.ya-distlimbs-0.05,tibiaR.ya-distlimbs-0.12-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHR1=data.frame(x=fingerHR1.x,y=fingerHR1.y)
  fingerHR2.x <- c(metatR.xc,metatR.xc+0.15+0.05*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHR2.y <- c(tibiaR.ya-distlimbs-0.025,tibiaR.ya-distlimbs-0.05-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHR2=data.frame(x=fingerHR2.x,y=fingerHR2.y)
  fingerHR3.x <- c(metatR.xc-0.1,metatR.xc+0.05*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHR3.y <- c(tibiaR.ya-distlimbs-0.025,tibiaR.ya-distlimbs-0.15-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2]))
  fingerHR3=data.frame(x=fingerHR3.x,y=fingerHR3.y)
  fingerFL1.x <- c(radiusL.xa,radiusL.xa)
  fingerFL1.y <- c(metacL.yb,metacL.yb+0.1+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFL2.x <- c(radiusL.xd,radiusL.xd)
  fingerFL2.y <- c(metacL.yb,metacL.yb+0.1+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFL3.x <- c((radiusL.xa+radiusL.xd)/2,(radiusL.xa+radiusL.xd)/2)
  fingerFL3.y <- c(metacL.yb,metacL.yb+0.1+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFL4.x <- c(radiusL.xa,radiusL.xa-0.05-0.05*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFL4.y <- c(radiusL.yb+0.02,radiusL.yb+0.01+0.07+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFL5.x <- c(radiusL.xd,radiusL.xd+0.05+0.05*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFL5.y <- c(radiusL.yb+0.02,radiusL.yb+0.01+0.07+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFR1.x <- c(radiusR.xa,radiusR.xa)
  fingerFR1.y <- c(metacL.yb,metacL.yb+0.1+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFR2.x <- c(radiusR.xd,radiusR.xd)
  fingerFR2.y <- c(metacL.yb,metacL.yb+0.1+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFR3.x <- c((radiusR.xa+radiusR.xd)/2,(radiusR.xa+radiusR.xd)/2)
  fingerFR3.y <- c(metacL.yb,metacL.yb+0.1+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFR4.x <- c(radiusR.xa,radiusR.xa-0.05-0.05*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFR4.y <- c(radiusL.yb+0.02,radiusL.yb+0.01+0.07+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFR5.x <- c(radiusR.xd,radiusR.xd+0.05+0.05*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  fingerFR5.y <- c(radiusL.yb+0.02,radiusL.yb+0.01+0.07+0.5*(X*prueba.pca$rotation[14]+Y*prueba.pca$rotation[14,2]))
  #fingerFL4=data.frame(x=fingerFL4.x,y=fingerFL4.y)
  #fingerFL5=data.frame(x=fingerFL5.x,y=fingerFL5.y)
  b <- ggplot(df) + geom_point() + xlim(-1.2, 1.2) + ylim(-1, 1) +
    geom_polygon(data=head, mapping=aes(x=x, y=y)) +
    geom_polygon(data=body, mapping=aes(x=x, y=y)) +
    geom_polygon(data=humerusL, mapping=aes(x=x, y=y)) +
    geom_polygon(data=radiusL, mapping=aes(x=x, y=y)) +
    geom_polygon(data=humerusR, mapping=aes(x=x, y=y)) +
    geom_polygon(data=radiusR, mapping=aes(x=x, y=y)) +
    geom_polygon(data=femurL, mapping=aes(x=x, y=y)) +
    geom_polygon(data=femurR, mapping=aes(x=x, y=y)) +
    geom_polygon(data=tibiaL, mapping=aes(x=x, y=y)) +
    geom_polygon(data=tibiaR, mapping=aes(x=x, y=y)) +
    geom_polygon(data=metacL, mapping=aes(x=x, y=y)) +
    geom_polygon(data=metacR, mapping=aes(x=x, y=y)) +
    geom_polygon(data=metatL, mapping=aes(x=x, y=y)) +
    geom_polygon(data=metatR, mapping=aes(x=x, y=y)) +
    geom_segment(aes(x=fingerHL1.x[1], y=fingerHL1.y[1], xend=fingerHL1.x[2], yend=fingerHL1.y[2]),size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHL1.x[2], y=fingerHL1.y[2], xend=fingerHL1.x[2], yend=fingerHL1.y[2]-0.1-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2])), size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHL2.x[1], y=fingerHL2.y[1], xend=fingerHL2.x[2], yend=fingerHL2.y[2]),size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHL2.x[2], y=fingerHL2.y[2], xend=fingerHL2.x[2], yend=fingerHL2.y[2]-0.1-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2])), size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHL3.x[1], y=fingerHL3.y[1], xend=fingerHL3.x[2], yend=fingerHL3.y[2]),size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHL3.x[2], y=fingerHL3.y[2], xend=fingerHL3.x[2], yend=fingerHL3.y[2]-0.1-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2])), size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHR1.x[1], y=fingerHR1.y[1], xend=fingerHR1.x[2], yend=fingerHR1.y[2]),size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHR1.x[2], y=fingerHR1.y[2], xend=fingerHR1.x[2], yend=fingerHR1.y[2]-0.1-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2])), size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHR2.x[1], y=fingerHR2.y[1], xend=fingerHR2.x[2], yend=fingerHR2.y[2]),size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHR2.x[2], y=fingerHR2.y[2], xend=fingerHR2.x[2], yend=fingerHR2.y[2]-0.1-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2])), size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHR3.x[1], y=fingerHR3.y[1], xend=fingerHR3.x[2], yend=fingerHR3.y[2]),size=1.1,col="grey20") +
    geom_segment(aes(x=fingerHR3.x[2], y=fingerHR3.y[2], xend=fingerHR3.x[2], yend=fingerHR3.y[2]-0.1-0.5*(X*prueba.pca$rotation[10]+Y*prueba.pca$rotation[10,2])), size=1.1,col="grey20") +
    geom_segment(aes(x=fingerFL1.x[1], y=fingerFL1.y[1], xend=fingerFL1.x[2], yend=fingerFL1.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFL2.x[1], y=fingerFL2.y[1], xend=fingerFL2.x[2], yend=fingerFL2.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFL3.x[1], y=fingerFL3.y[1], xend=fingerFL3.x[2], yend=fingerFL3.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFL4.x[1], y=fingerFL4.y[1], xend=fingerFL4.x[2], yend=fingerFL4.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFL5.x[1], y=fingerFL5.y[1], xend=fingerFL5.x[2], yend=fingerFL5.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFR1.x[1], y=fingerFR1.y[1], xend=fingerFR1.x[2], yend=fingerFR1.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFR2.x[1], y=fingerFR2.y[1], xend=fingerFR2.x[2], yend=fingerFR2.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFR3.x[1], y=fingerFR3.y[1], xend=fingerFR3.x[2], yend=fingerFR3.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFR4.x[1], y=fingerFR4.y[1], xend=fingerFR4.x[2], yend=fingerFR4.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFR5.x[1], y=fingerFR5.y[1], xend=fingerFR5.x[2], yend=fingerFR5.y[2]), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFL4.x[2], y=fingerFL4.y[2], xend=fingerFL4.x[2], yend=fingerFL1.y[2]-0.04), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFL5.x[2], y=fingerFL5.y[2], xend=fingerFL5.x[2], yend=fingerFL1.y[2]-0.04), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFR4.x[2], y=fingerFR4.y[2], xend=fingerFR4.x[2], yend=fingerFR1.y[2]-0.04), size=1.05,col="grey20") +
    geom_segment(aes(x=fingerFR5.x[2], y=fingerFR5.y[2], xend=fingerFR5.x[2], yend=fingerFR1.y[2]-0.04), size=1.05,col="grey20") +
    geom_curve(aes(x = body.xe-0.01, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xe-0.03, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xe-0.05, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xe-0.07, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xe-0.09, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xe-0.11, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xe-0.13, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xe-0.15, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xe-0.17, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.01, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.03, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.05, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.07, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.09, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.11, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.13, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.15, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.17, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = body.xf+0.17, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = 0, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = 0.02, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = -0.02, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = 0.04, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_curve(aes(x = -0.04, y = -0.53, xend = 1, yend = -1), size=2, col="grey20") +
    geom_segment(aes(x=body.xe, y=-0.53, xend=body.xf, yend=-0.53), size=2,col="grey20") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) + theme_void() 
  
  c <- ggarrange(b,a,ncol = 2, nrow = 1, widths=c(1,1.4))
  print(c)
  dev.off()}


png_files <- list.files("~/Desktop/R/", pattern = ".*png$", full.names = TRUE)
png_files <- mixedsort(png_files, decreasing=FALSE)
gifski(png_files, gif_file = "animation.gif", width = 800, height = 400, delay = 0.1)

#file.remove(list.files(pattern=".png"))
