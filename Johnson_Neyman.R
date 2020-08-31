setwd("/Users/Ken/Desktop/R")
data <- read.table("rawdata_thoracicus.txt",header=T, na.strings = c("x", "NA"))
data3 <- read.table("microlophus2.txt",header=T, na.strings = c("x", "NA"))

data3 <- data[- which(data$age=="j"), ]
data3 <- data3[- which(data3$species2=="thor_j"), ]
model <- aov(log(data3$hh)~log(data3$hl)*data3$species2)
TukeyHSD(model,"data3$species2")

summary(aov(log(data3$hh)~log(data3$hl)*data3$species2))
data <- data3
tiff(filename = "/Users/KenToyama/Desktop/without.tiff", 
     width = 1600, height = 1300, units= "px", res = 300)

dev.off() 

plot((data3$svl),(data3$head_height), ylab="Head width", xlab='svl')
polygon(c(xlower,xlower,xupper,xupper),c(0,3.2,3.2,0),
        col=rgb(224, 224, 224,maxColorValue=255,alpha=130), border=NA)
abline(v=xlower,lty=1, lwd=0.5)
abline(v=xupper,lty=1, lwd=0.5)
dat <- data[ which(data$sex=="0"), ]
points(log(dat$head_length),log(dat$head_height), pch=16 ,col="black")
newdata <- data[ which(data$sex=='1' ), ]
#clip(4.0,4.5, -100, 100)
abline(lm(log(newdata$head_height) ~ log(newdata$head_length)))
newdata <- data[ which(data$sex=='0' ), ]
#points(log(newdata$hl),log(newdata$hh), pch=16, col="white")
#clip(3.4,4.0, -100, 100)
abline(lm(log(newdata$head_height) ~ log(newdata$head_length)), lty=2)

group1 <- data[ which(data$sex=="1"), ]
group2 <- data[ which(data$sex=="0"), ]

reg1 <- lm(log(adults$head_width)~log(adults$svl))
reg2 <- lm(log(juveniles$head_width)~log(juveniles$svl))
coef1 <- reg1$coefficients
coef1[2]
###########################################
n1 <- 26
n2 <- 29
X1 <- group1$head_length
X2 <- group2$head_length
Y1 <- group1$head_height
Y2 <- group2$head_height
reg1 <- lm(log(Y1)~log(X1))
reg2 <- lm(log(Y2)~log(X2))
coef1 <- reg1$coefficients
coef2 <- reg2$coefficients
Fcvalue <- qf(.95, df1=1, df2=n1+n2-4)
xmean1 <- mean(log(X1))
xmean2 <- mean(log(X2))
xmeansq1 <- mean((log(X1))^2)
xmeansq2 <- mean((log(X2))^2)
sumx1 <- 0
sumx2 <- 0
sumy1 <- 0
sumy2 <- 0
sumxy1 <- 0
sumxy2 <- 0

xcoord_group1 <- log(X1)    # X1
ycoord_group1 <- log(Y1)    # Y1

xcoord_group2 <- log(X2)    # X2
ycoord_group2 <- log(Y2)    # Y2

zx1 <- ((sum(xcoord_group1))^2)/n1
zx2 <- ((sum(xcoord_group2))^2)/n2
zy1 <- ((sum(ycoord_group1))^2)/n1
zy2 <- ((sum(ycoord_group2))^2)/n2
zxy1 <- ((sum(xcoord_group1))*(sum(ycoord_group1)))/n1
zxy2 <- ((sum(xcoord_group2))*(sum(ycoord_group2)))/n2

######## sumx1 ########
c <- 0
while (c<n1) {
  sumx1 <- sumx1 + (((xcoord_group1[c+1])^2) - zx1)
  c <- c+1
}
sumx1

######## sumx2 ########
c <- 0
while (c<n2) {
  sumx2 <- sumx2 + (((xcoord_group2[c+1])^2) - zx2)
  c <- c+1
}
sumx2

######## sumy1 ########
c <- 0
while (c<n1) {
  sumy1 <- sumy1 + (((ycoord_group1[c+1])^2) - zy1)
  c <- c+1
}
sumy1

######## sumy2 ########
c <- 0
while (c<n2) {
  sumy2 <- sumy2 + (((ycoord_group2[c+1])^2) - zy2)
  c <- c+1
}
sumy2

######## sumxy1 ########
c <- 0
while (c<n1) {
  sumxy1 <- sumxy1 + (((xcoord_group1[c+1])*(ycoord_group1[c+1]))-zxy1)
  c <- c+1
}
sumxy1

######## sumxy2 ########
c <- 0
while (c<n2) {
  sumxy2 <- sumxy2 + (((xcoord_group2[c+1])*(ycoord_group2[c+1]))-zxy2)
  c <- c+1
}
sumxy2

########################
########################
########################
######## SSres ########

SSres <- (sumy1-(((sumxy1)^2)/sumx1))+(sumy2-(((sumxy2)^2)/sumx2))
SSres

########################
########################
########################
######## ABC ########
a1 <- coef1[1]
a2 <- coef2[1]
b1 <- coef1[2]
b2 <- coef2[2]
A <- (-Fcvalue/(n1+n2-4))*(SSres)*((1/sumx1)+(1/sumx2))+((b1-b2)^2)
B <- (Fcvalue/(n1+n2-4))*(SSres)*((xmean1/sumx1)+(xmean2/sumx2))+((a1-a2)*(b1-b2))
C <- (-Fcvalue/(n1+n2-4))*(SSres)*(((n1+n2)/(n1*n2))+(xmeansq1/sumx1)+(xmeansq2/sumx2))+((a1-a2)^2)

########################

xlower <- (-B-sqrt((B^2)-A*C))/A
xupper <- (-B+sqrt((B^2)-A*C))/A

xlower
xupper

