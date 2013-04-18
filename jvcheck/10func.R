library(MASS)
library(foreach)
gendata <- function(n,m,theta,ar,x.innov,e.innov) {
    vv <- x.innov((n+200)*m)    
    xx <- ts(filter(vv,filter=ar,method="recursive",sides=1),frequency=m)   
   
    ms <- msim(n,theta,xx,e.innov)
    x <- window(xx,start=start(ms$y))
    list(y=ms$y,x=x,vv=vv,e=ms$e)
}

msim <- function(n,theta,x,e.innov) {
    m <- frequency(x)
    n.x <- length(x)
    
    X <- fmls(x,length(theta)-1,m)
    xt <- as.vector(X%*%theta)
    
    xet <- xt+e.innov(length(xt))
    list(y=ts(xet[nrow(X)-n:1+1],start=end(x)[1]-n+1,frequency=1),e=xet-xt)    
}

genhahr <- function(N,n,hh,ar,inno.x,inno.e) {
    require(foreach)
    rr <- foreach(i=1:N,.combine="rbind") %do% {
        aa <- gendata(n,12,nealmon(hh,12),ar,inno.x,inno.e)       
        bb <- midas_r(y~fmls(x,11,12,nealmon),data.frame(y=aa$y),data.frame(x=aa$x),start=list(x=hh+rnorm(length(hh))/10))
        PHI=vcovHAC(bb$unrestricted,sandwich=FALSE)
        res1 <- hAhr.test(bb,PHI=PHI)
        res2 <- hAhrfix.test(bb,PHI=PHI)
        c(res1$statistic,res2$statistic)
    }
    rr
}


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
