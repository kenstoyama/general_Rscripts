install.packages("plotrix")
install.packages("grDevices")
install.packages("gtools")
install.packages("gifski")
library(plotrix)
library(grDevices)
library(gtools)
library(gifski)
setwd("~/Desktop/GitHub/general_Rscripts/")
data <- read.csv("avg_gif.csv")
data <- data[order(-data$femur),]
gray <- gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = FALSE)
color <- c(rep(gray[1],9),rep(gray[2],9),rep(gray[3],9),rep(gray[4],9),rep(gray[5],9),
  rep(gray[6],9),rep(gray[7],9),rep(gray[8],9),rep(gray[9],9),rep(gray[10],9))
data <- cbind(data,color)
c <- 1

while (c<91) {
  png(paste("test",c,".png"), units="in", width=6, height=5, res=300)
  y <- data[c,4]
  x <- data[c,3]
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot(0:5,xlim=c(-2,2),ylim=c(-4,4),type="n",xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")
  draw.circle(0,0,radius = 0.75,nv=100,col=y,border = "black")
  draw.circle(0,0,radius = 0.75*x,nv=100,col="white",border = "black")
  #xleft, ybottom, xright, ytop
  #rect(-1,2.5,1,3,col="black")
  rect(-1,2.5,-0.8,3,col=gray[10],lwd=0)
  rect(-0.8,2.5,-0.6,3,col=gray[9],lwd=0)
  rect(-0.6,2.5,-0.4,3,col=gray[8],lwd=0)
  rect(-0.4,2.5,-0.2,3,col=gray[7],lwd=0)
  rect(-0.2,2.5,0,3,col=gray[6],lwd=0)
  rect(0,2.5,0.2,3,col=gray[5],lwd=0)
  rect(0.2,2.5,0.4,3,col=gray[4],lwd=0)
  rect(0.4,2.5,0.6,3,col=gray[3],lwd=0)
  rect(0.6,2.5,0.8,3,col=gray[2],lwd=0)
  rect(0.8,2.5,1,3,col=gray[1],lwd=0)
  text(-1.2,2.8,"low")
  text(1.25,2.8,"high")
  text(0,3.5,"bone mineral density")
  spp_name <- paste("Anolis",data[c,1])
  text(0,-2.5,spp_name,font=4)
  c <- c+1
  print(c)
  dev.off()
}
png_files <- list.files("~/Desktop/GitHub/general_Rscripts/", pattern = ".*png$", full.names = TRUE)
png_files <- mixedsort(png_files, decreasing=FALSE)
gifski(png_files, gif_file = "BMD_gif.gif", width = 600, height = 500, delay = 0.2)


file.remove(list.files(pattern=".png"))


