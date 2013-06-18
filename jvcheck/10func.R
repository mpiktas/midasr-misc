library(MASS)
library(foreach)
library(sandwich)

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

gencr <- function(N,n,wfun,wpar,wfun1,wpar1=wpar,dk,inno.x,inno.e,m=12,ar=0.6,...) {
    require(foreach)
    fr <- formula(y~fmls(x,a,b,c))
    fr[[c(3,3)]] <- dk-1
    fr[[c(3,4)]] <- m
    fr[[c(3,5)]] <- as.name(wfun1)
    
    rr <- foreach(i=1:N,.combine="c",.errorhandling="pass") %do% {
        aa <- gendata(n,m,wfun(wpar,dk),ar,inno.x,inno.e)        
        bb <- try(midas_r(fr,data.frame(y=aa$y),data.frame(x=aa$x),start=list(x=wpar1+rnorm(length(wpar1))/10)))
        ac <- function() {
            for(badi in 1:3) {
     #           print(badi)
                st <- wpar1+rnorm(length(wpar1))/5             
                bb <- try(midas_r(fr,data.frame(y=aa$y),data.frame(x=aa$x),start=list(x=st)))
                if(!inherits(bb,"try-error")) {
                    conv <- bb$opt$convergence
                    if(conv==0) break                    
                }               
            }
            bb
        }
        if(inherits(bb,"try-error")){          
            bb <- ac()            
        }
        else {            
            if(bb$opt$convergence!=0) bb <- ac()            
        }
        if(inherits(bb,"try-error")) {
            list(bb)
        }
        else {
            PHI <- vcovHAC(bb$unrestricted,sandwich=FALSE)
            res1 <- hAhr.test(bb,PHI)
            res2 <- hAhrfix.test(bb,PHI)
            list(list(orig=res1,fix=res2))
        }
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

fixrr1 <- function(x) {
    
}
getres4 <- function(rr1,what=c("orig","fix")) {
    what <- match.arg(what)
    require(data.table)
    require(reshape2)
    uu <- lapply(rr1,function(y){
        y$gen$message <- NULL
        c(y$l$i,y$l$n,y$l$N,y$l$dk,
          sapply(y$gen,function(x) {
              if(is.null(names(x))) {
                  NA
              } else {
                  x[[what]]$statistic
              }
          }))        
    })
    ee <- melt(data.frame(do.call("rbind",uu)),id=1:4)
    ee <- ee[,-c(3,5)]
    ee <- data.table(ee)
    setnames(ee,c("V1","V2","V4"), c("i","n","dk"))
    setkey(ee,i,n,dk)
    ee
}

calcp <- function(ee,thres=7.81) {
    r1 <- ee[,sum(na.omit(value>thres)),by=key(ee)]
    r2 <- ee[,sum(!is.na(value)),by=key(ee)]
    r3 <- ee[,sum(na.omit(value<0)),by=key(ee)]
    fr <- r1[r2[r3]]
    fr$p <- fr$V1/fr$V1.2
    fr
}
