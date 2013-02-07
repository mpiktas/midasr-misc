require(midasr)

imidas_r_fast <- midasr:::imidas_r_fast
  
gendata <- function(n,dk,m,ar,innov.sd,theta) {
   # vv <- simplearma.sim(model=list(ar=ar),(n+10)*m,innov.sd,m)
    vv <- ts(arima.sim((n+10)*m,model=list(ar=ar),sd=1),frequency=m)
    xx <- ts(cumsum(vv),start=start(vv),frequency=frequency(vv))
    y <- midas.sim(n,theta,xx,1)
    x <- window(xx,start=start(y))
    list(y=y,x=x,vv=vv)
}

mm.sim <- function(n,dk,m,ar,innov.sd,theta) {
    vv <- ts(arima.sim((n+10)*m,model=list(ar=ar),sd=1),frequency=m)
    xx <- ts(cumsum(vv),start=start(vv),frequency=frequency(vv))
    X <- fmls(xx,length(theta)-1,m)
    xt <- as.vector(X%*%theta)
    e <- rnorm(n,sd=innov.sd)
    y<- xt[nrow(X)-n:1+1]+e
    list(vv=vv,xx=xx,X=X,xt=xt,e=e,y=y)
}

size.imr <- function(n,dk,m,ar,innov.sd=1,weight,cf,simplify=TRUE) {
    theta <- weight(cf,dk+1)
    dt <- gendata(n,dk,m,ar,innov.sd,theta)
    y <- dt$y
    x <- dt$x
    imr.nls <- imidas_r(y~fmls(x,dk,m,weight)-1,model="twosteps",start=list(x=rnorm(length(cf))),control=list(maxit=1000,reltol=sqrt(.Machine$double.eps)/10,hessian=T))
    imr.td <-  imidas_r(y~fmls(x,dk,m,weight)-1,model="reduced",start=list(x=rnorm(length(cf))),control=list(maxit=1000,reltol=sqrt(.Machine$double.eps)/10,hessian=T))
    
    if(simplify) {
        inls <- ihAh.test(imr.nls)
        itd <- ihAh.test(imr.td)
        nm <- c("statistic","p.value")
        return(c(unlist(inls[nm]),unlist(itd[nm])))
    }
    else {
        list(nls=imr.nls,td=imr.td,theta=theta,data=dt)
    }
}

size.imr.fast <- function(n,dk,m,ar,innov.sd=1,weight,cf,simplify=TRUE) {
    theta <- weight(cf,dk+1)
    dt <- gendata(n,dk,m,ar,innov.sd,theta)
    y <- dt$y
    x <- dt$x

    V <- dmls(x,dk,m)
    xmd <- mls(x,c(dk,dk+1),m)    
    y <- as.numeric(y)
    model <- na.omit(cbind(y,xmd[,2],V))
    model1 <- na.omit(cbind(y,xmd[,1],V[,1:dk]))
    
    imr.nls <- imidas_r_fast(y,x,dk,weight,start=rnorm(length(cf)),imodel="twosteps",model.matrix=model)
    imr.td <- imidas_r_fast(y,x,dk,weight,start=rnorm(length(cf)),imodel="reduced",model.matrix=model1)
       
    if(simplify) {
        inls <- ihAh.nls.test(imr.nls)
        itd <- ihAh.Td.test(imr.td)
        nm <- c("statistic","p.value")
        return(c(unlist(inls[nm]),unlist(itd[nm])))
    }
    else {
        list(nls=imr.nls,td=imr.td,theta=theta,data=dt,model=model,model1=model1)
    }
}

sim.rowdata <- function(N,n,dk,m,ar,innov.sd=1,weight,cf) {
    res <- foreach(i=1:N,.combine="c",.errorhandling="pass") %dopar% {
        list(try(size.imr.fast(n,dk,m,ar,weight=weight,cf=cf,simplify=FALSE)))
    }
}

sim.rowdata2 <- function(N,n,dk,m,ar,innov.sd=1,weight,cf,simplify=FALSE) {
    res <- vector("list",N)
    for(i in 1:N){
        res[[i]] <- try(size.imr.fast(n,dk,m,ar,weight=weight,cf=cf,simplify=simplify))
    }
    res
}
