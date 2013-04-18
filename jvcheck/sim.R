library(midasr)
library(foreach)
source("10func.R")
source("hahr.R")
##The parameter function

h01 <- c(10,0,0)
h02 <- c(10,-5,-5)
h03 <- c(10,1,1)
h04 <- c(10,-0.5,0.04)
h05 <- c(10,0.5,-0.04)
h06 <- c(10,0.005,0.02)
##Constant to multiply to adjust to our nealmon function
const <- c(1,100,100^2)
HH<-rbind(h01,h02,h03,h04,h05,h06)


xino <- function(n)rnorm(n)


rr1<-foreach(i=4:6,.combine="c",.errorhandling="pass") %dopar% {
   ii<-foreach(n=c(300,1000,2000),.combine="c",.errorhandling="pass") %do% {
     list(genhahr(1,n,HH[i,]*const,ar=0.6,xino,xino))
   }
   list(ii)
}

if(FALSE) {
n <- 500
vals <- matrix(NA, nrow = n, ncol = 3)
freqs <- c(300, 1000, 2000)

for(d in 1:3){
  for(i in 1:nrow(vals)){
    x <- simplearma.sim(list(ar=0.6),2500 * 12,1,12)
    y <- midas.sim(freqs[d],theta0,x,1)
    x <- window(x,start=start(y))
    mr <- midas_r(y~fmls(x,4*12-1,12,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(x=c(-0.1,0.1,-0.1,-0.001)),dk=4*12)
    vals[i, d] <- hAhrfix.test(mr)$statistic
  }
}
save(vals,file="valsfix.RData")
}


#plot(seq(0, 100, 0.1), 
#     dchisq(seq(0, 100, 0.1), 44), 
#     xlab = "red - 300 obs., magenta - 1000, blue - 2000", 
#3     ylab= "")
#lines(density(vals[,1]), col = "red")
#lines(density(vals[,2]), col = "magenta")
#lines(density(vals[,3]), col = "blue")
