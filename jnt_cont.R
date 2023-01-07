data <- avg_M
X <- "svl_male"
Y <- "SSD"
g <- "femur"

jnt_cont <- function(X,Y,g,data,phylo=F,tree){
  mod <- summary(lm(data[,Y]~data[,X]*data[,g]))
  mod.out <- mod
  if(phylo==T){
    Xi <- data[,X]
    Yi <- data[,Y]
    gi <- data[,g]
    mod <- gls(Yi~Xi*gi, correlation=corPagel(1,tree),
               na.action = na.omit)
    mod.out <- summary(mod)
  }
  varcov <- vcov(mod.out)
  coef1 <- mod.out$coefficients[2]
  coef3 <- mod.out$coefficients[4]
  var_coef1 <- varcov[2,2]
  var_coef3 <- varcov[4,4]
  cov_coefs1_3 <- varcov[2,4]
  
  a <- (1.96^2)*(var_coef3)-(coef3)^2
  b <- 2*((1.96^2)*(cov_coefs1_3)-(coef1*coef3))
  c <- (1.96^2)*var_coef1-(coef1)^2
  
  x1 <- (-b-sqrt((b^2)-(4*a*c)))/(2*a)
  x2 <- (-b+sqrt((b^2)-(4*a*c)))/(2*a)
  
  m <- c("It was not possible to calculate regions of non-significance. The difference between slopes might not be statistically significant")
  if(((b^2)-4*a*c)<0){
    warning(m)
  }
  
  plot(data[,X],data[,Y],xlab=X,ylab=Y)
  
  val <- min(min(data[,g],na.rm = T),min(x1,x2))
  valmin <- min(data[,g],na.rm = T)
  val_max <- max(max(data[,g],na.rm = T),max(x1,x2))
  valmax <- max(data[,g],na.rm = T)
  nn <- (val_max-val)/100
  aaa <- c()
  c <- 1
  while(val+nn<val_max){
    if(val+nn>x2 & val+nn<x1){
      abline(a=(mod.out$coefficients[1]+mod.out$coefficients[3]*val),b=(mod.out$coefficients[2]+mod.out$coefficients[4]*val),col=alpha("grey",0.5))
      aaa[c] <- val
      val <- val+nn
      c <- c+1
    } else {
      abline(a=(mod.out$coefficients[1]+mod.out$coefficients[3]*val),b=(mod.out$coefficients[2]+mod.out$coefficients[4]*val),col=alpha("lightblue",0.5))
      aaa[c] <- val
      val <- val+nn
      c <- c+1
    }
  }
  abline(a=(mod.out$coefficients[1]+mod.out$coefficients[3]*min(data[,g],na.rm = T)),
         b=(mod.out$coefficients[2]+mod.out$coefficients[4]*min(data[,g],na.rm = T)),
         col="black",lwd=1,lty=2)
  abline(a=(mod.out$coefficients[1]+mod.out$coefficients[3]*max(data[,g],na.rm = T)),
         b=(mod.out$coefficients[2]+mod.out$coefficients[4]*max(data[,g],na.rm = T)),
         col="black",lwd=1,lty=2)
  results <- list("coeff" = mod.out,"lower limit of moderator" = x1,
                  "upper limit of moderator" = x2,
                  "lower data limit" = valmin, "upper data limit" = valmax)
  return(results)
}

jnt_cont("svl_male","femur","SSD",data=avg_M,phylo=F,tree=phy_M)
jnt_cont("svl_male","final_rad","SSD",data=avg_M,phylo=F,tree=phy_M)
jnt_cont("svl","head_height","head_length",data=thor,phylo=F,tree=phy_M)
mod <- lm(femur~svl_male*SSD, data=avg_M)
mod <- lm(head_height~svl*head_length, data=thor)
johnson_neyman(model=mod,pred=svl_male,modx=SSD)

thor$svl <- log(thor$svl)
thor$head_length <- log(thor$head_length)
thor$head_height <- log(thor$head_height)
mod.out
mod <- summary(lm(thor$head_height~thor$head_length*thor$svl))

