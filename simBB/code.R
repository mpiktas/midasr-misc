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

size.imr.fast <- function(n,dk,m,ar,innov.sd=1,weight0,cf0,weight1=weight0,cf1=cf0,simplify=TRUE) {
    theta <- weight0(cf0,dk+1)
    dt <- gendata(n,dk,m,ar,innov.sd,theta)
    y <- dt$y
    x <- dt$x

    V <- dmls(x,dk,m)
    xmd <- mls(x,c(dk,dk+1),m)    
    y <- as.numeric(y)
    model <- na.omit(cbind(y,xmd[,2],V))
    model1 <- na.omit(cbind(y,xmd[,1],V[,1:dk]))
    
    imr.nls <- imidas_r_fast(y,x,dk,weight1,start=cf1+rnorm(length(cf1))/5,imodel="twosteps",model.matrix=model)
    imr.td <- imidas_r_fast(y,x,dk,weight1,start=cf1+rnorm(length(cf1))/5,imodel="reduced",model.matrix=model1)
    imr.td2 <- imidas_r_fast(y,x,dk,weight1,start=cf1+rnorm(length(cf1))/5,imodel="reduced2",model.matrix=model1)
    
    if(simplify) {
        inls <- ihAh.nls.test(imr.nls)
        itd <- ihAh.Td.test(imr.td)
        itd2 <- ihAh.Td.test(imr.td2)
        nm <- c("statistic","p.value")
        return(list(t.nls=inls[nm],t.td=itd[nm],t.td2=itd2[nm],nls=imr.nls$opt,td=imr.td$opt,td2=imr.td2$opt))
    }
    else {
        list(nls=imr.nls,td=imr.td,td2=imr.td2,theta=theta,data=dt,model=model,model1=model1)
    }
}

sim.rowdata <- function(N,n,dk,m,ar,innov.sd=1,weight,cf) {
    res <- foreach(i=1:N,.combine="c",.errorhandling="pass") %dopar% {
        list(try(size.imr.fast(n,dk,m,ar,weight=weight,cf=cf,simplify=FALSE)))
    }
}

sim.rowdata2 <- function(N,n,dk,m,ar,innov.sd=1,weight0,cf0,weight1=weight0,cf1=cf0,simplify=FALSE) {
    res <- vector("list",N)
    for(i in 1:N){
        res[[i]] <- try(size.imr.fast(n,dk,m,ar,weight0=weight0,cf0=cf0,weight1=weight1,cf1=cf1,simplify=simplify))
    }
    res
}

f.theta.kz3 <- function(p, dk) {
    i <- (1:dk)/100
    (p[1]*i)*exp(p[2]*i + p[3]*i^2)
}

simtb <- function(N,param,weight0,cf0,weight1=weight0,cf1=cf0,innov.sd=1,simplify=FALSE) {    
    foreach(i=1:nrow(param),.combine="c",.errorhandling="pass") %dorng% {
        n <- param$n[i]
        dk <- param$d[i]
        m <- param$m[i]
        ar <- param$ar[i]
        list(try(sim.rowdata2(N,n,dk,m,ar,innov.sd=innov.sd,weight0,cf0,weight1,cf1,simplify=TRUE)))
    }
}

selstat <- function(x,type="td",thresh=0.1,stat="p.value") {
    st <- paste("t",type,sep=".")
    if(is.list(x)) {
	test <- sqrt(sum(x[[type]][["gr0"]]^2))
	if(test<thresh& x[[type]]$convergence==0) return(x[[st]][[stat]])
    }
    NA
}
calccol <- function(tb,type="td",thresh=0.1) {
    lc <- sapply(tb,function(l)sapply(l,selstat,type=type,thresh=thresh))
    res<-apply(lc,2,function(x){
	xx<-na.omit(x)
	c(sum(xx<0.05)/length(xx),length(xx))
    })
    res<-t(res)
    colnames(res) <- c(type,paste("n",type,sep=""))
    res
}
cstarhat<-function(row,alpha=0.05,type="td",thresh=0.1) {
    t0<-sapply(row,selstat,type=type,thresh=thresh,stat="statistic")
    t0 <- na.omit(t0)
    quantile(t0,1-alpha)
}
adjpow <- function(tb,pow,alpha=0.05,type="td",thresh=0.1) {
    cstar<-sapply(tb,cstarhat,alpha=alpha,type=type,thresh=0.1)
    lc <- lapply(pow,function(l)sapply(l,selstat,type=type,thresh=thresh,stat="statistic"))
    res <- mapply(function(row,cs){
	xx<- na.omit(row)
	c(sum(xx>cs)/length(xx),length(xx))
    },lc,cstar)
    res <- t(res) 
    colnames(res) <- c(type,paste("n",type,sep=""))
    res
}

sizetable<-function(tbd,param,thresh=0.1) {
    snls <- calccol(tbd,type="nls",thresh=thresh)
    std <- calccol(tbd,type="td",thresh=thresh)
    std2 <- calccol(tbd,type="td2",thresh=thresh)
    cbind(param,snls,std,std2)
}

aptable<-function(tbd,pow,param,alpha=0.05,thresh=0.1) {
    anls<-adjpow(tbd,pow,alpha,type="nls",thresh=thresh)
    astd <- adjpow(tbd,pow,alpha,type="td",thresh=thresh)
    astd2<-adjpow(tbd,pow,alpha,type="td2",thresh=thresh)
    cbind(param,anls,astd,astd2)
}
