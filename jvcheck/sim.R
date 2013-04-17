library(midasr)

##The parameter function
theta.h0 <- function(p, dk) {
  i <- (1:dk-1)
  (p[1] + p[2]*i)*exp(p[3]*i + p[4]*i^2)
}

##Generate coefficients
theta0 <- theta.h0(c(-0.1,0.1,-0.1,-0.001),4*12)

##The gradient function
grad.h0<-function(p, dk) {
  i <- (1:dk-1)
  a <- exp(p[3]*i + p[4]*i^2)
  cbind(a, a*i, a*i*(p[1]+p[2]*i), a*i^2*(p[1]+p[2]*i))
}

n <- 500
vals <- matrix(NA, nrow = n, ncol = 3)
freqs <- c(300, 1000, 2000)
for(d in 1:3){
  for(i in 1:nrow(vals)){
    x <- simplearma.sim(list(ar=0.6),2500 * 12,1,12)
    y <- midas.sim(freqs[d],theta0,x,1)
    x <- window(x,start=start(y))
    mr <- midas_r(y~fmls(x,4*12-1,12,theta.h0)-1,data.frame(y=y),data.frame(x=x),start=list(x=c(-0.1,0.1,-0.1,-0.001)),dk=4*12)
    vals[i, d] <- hAhr.test(mr)$statistic
  }
}

plot(seq(0, 100, 0.1), 
     dchisq(seq(0, 100, 0.1), 44), 
     xlab = "red - 300 obs., magenta - 1000, blue - 2000", 
     ylab= "")
lines(density(vals[,1]), col = "red")
lines(density(vals[,2]), col = "magenta")
lines(density(vals[,3]), col = "blue")
