h01 <- c(10,0,0)
h02 <- c(10,-5,-5)
h03 <- c(10,1,1)
h04 <- c(10,-0.5,0.04)
h05 <- c(10,0.5,-0.04)
h06 <- c(10,0.005,0.02)
##Constant to multiply to adjust to our nealmon function
const <- c(1,100,100^2)
HH<-rbind(h01,h02,h03,h04,h05,h06)

theta.h0 <- function(p, dk) {
  i <- (1:dk-1)
  (p[1] + p[2]*i)*exp(p[3]*i + p[4]*i^2)
}

shfun <- list(nealmon1=nealmon,
              nealmon2=nealmon,
              nealmon3=nealmon,
              theta.h0=theta.h0
              )

shfname <- list(nealmon1="nealmon",
               nealmon2="nealmon",
               nealmon3="nealmon",
               theta.h0="theta.h0" 
              )

shpar <- list(nealmon1=HH[4,]*const,
              nealmon2=HH[5,]*const,
              nealmon3=HH[6,]*const,
              theta.h0=c(-0.1,0.1,-0.1,-0.001))

shdk <- list(nealmon1=12,
             nealmon2=12,
             nealmon3=12,
             theta.h0=48)


if(FALSE) {
    ##Plot the shapes
    library(ggplot2)
    library(reshape2)
    print(qplot(x=x,y=value,data=do.call("rbind",mapply(function(x,y,z,d)data.frame(x=1:z,value=x(y,z),name=d),shfun,shpar,shdk,as.list(names(shfun)),SIMPLIFY=FALSE)),geom="line")+facet_wrap(~name,scales="free"))
    
}
