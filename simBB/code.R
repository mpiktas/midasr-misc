require(midasr)

imidas_r_fast <- midasr:::imidas_r_fast
ihAh.Td.test <- midasr:::ihAh.Td.test
ihAh.nls.test <- midasr:::ihAh.nls.test
  
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
    imr.nls2 <- imidas_r_fast(y,x,dk,weight1,start=cf1+rnorm(length(cf1))/5,imodel="twosteps2",model.matrix=model)
    imr.td <- imidas_r_fast(y,x,dk,weight1,start=cf1+rnorm(length(cf1))/5,imodel="reduced",model.matrix=model1)
    imr.td2 <- imidas_r_fast(y,x,dk,weight1,start=cf1+rnorm(length(cf1))/5,imodel="reduced2",model.matrix=model1)
    
    if(simplify) {
        inls <- ihAh.nls.test(imr.nls,"ols")
        inls2 <- ihAh.nls.test(imr.nls2,"ols")
        ninls <- ihAh.nls.test(imr.nls,"nls")
        ninls2 <- ihAh.nls.test(imr.nls2,"nls")

        itd <- ihAh.Td.test(imr.td,"ols")
        itd2 <- ihAh.Td.test(imr.td2,"ols")
        nitd <- ihAh.Td.test(imr.td,"nls")
        nitd2 <- ihAh.Td.test(imr.td2,"nls")

        nm <- c("statistic","p.value")
        return(list(se.ols=list(
                        t.nls=inls[nm],
                        t.nls2=inls2[nm],
                        t.td=itd[nm],
                        t.td2=itd2[nm]),
                    se.nls=list(
                        t.nls=ninls[nm],
                        t.nls2=ninls2[nm],
                        t.td=nitd[nm],
                        t.td2=nitd2[nm]
                        ),
                    nls=imr.nls$opt,
                    nls2=imr.nls2$opt,
                    td=imr.td$opt,
                    td2=imr.td2$opt))
    }
    else {
        list(nls=imr.nls,nls2=imr.nls2,td=imr.td,td2=imr.td2,theta=theta,data=dt,model=model,model1=model1)
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
        if(is.matrix(cf1)) {
            cf1c <- cf1[,as.character(dk)]
        } else {
            cf1c <- cf1
        }
        list(try(sim.rowdata2(N,n,dk,m,ar,innov.sd=innov.sd,weight0,cf0,weight1,cf1c,simplify=TRUE)))
    }
}

selstat <- function(x,type="td",thresh=0.1,stat="p.value",se.type="ols") {
    st <- paste("t",type,sep=".")
    if(is.list(x)) {
        test <- sqrt(sum(x[[type]][["gr0"]]^2))
        if(test<thresh& x[[type]]$convergence==0)  {
            if(se.type %in% names(x)) {
                return(x[[se.type]][[st]][[stat]])
            } else {
                return(x[[st]][[stat]])
            }
        }
    }
    NA
}

calccol <- function(tb,type="td",thresh=0.1,se.type="ols") {
    lc <- sapply(tb,function(l)sapply(l,selstat,type=type,thresh=thresh,se.type=se.type))
    res<-apply(lc,2,function(x){
        xx<-na.omit(x)
        c(sum(xx<0.05)/length(xx),length(xx))
    })
    res<-t(res)
    colnames(res) <- c(type,paste("n",type,sep=""))
    res
}

cstarhat<-function(row,alpha=0.05,type="td",thresh=0.1,se.type="ols") {
    t0<-sapply(row,selstat,type=type,thresh=thresh,stat="statistic",se.type=se.type)
    t0 <- na.omit(t0)
    quantile(t0,1-alpha)
}

adjpow <- function(tb,pow,alpha=0.05,type="td",thresh=0.1,se.type="ols") {
    cstar<-sapply(tb,cstarhat,alpha=alpha,type=type,thresh=0.1,se.type=se.type)
    lc <- lapply(pow,function(l)sapply(l,selstat,type=type,thresh=thresh,stat="statistic",se.type=se.type))
    res <- mapply(function(row,cs){
	xx<- na.omit(row)
	c(sum(xx>cs)/length(xx),length(xx))
    },lc,cstar)
    res <- t(res) 
    colnames(res) <- c(type,paste("n",type,sep=""))
    res
}

sizetable<-function(tbd,param,thresh=0.1,se.type="se.ols") {
    std <- calccol(tbd,type="td",thresh=thresh,se.type=se.type)
    snls <- calccol(tbd,type="nls",thresh=thresh,se.type=se.type)
    snls2 <- calccol(tbd,type="nls2",thresh=thresh,se.type=se.type)
    std2 <- calccol(tbd,type="td2",thresh=thresh,se.type=se.type)
    cbind(param,snls,snls2,std,std2)
}

aptable<-function(tbd,pow,param,alpha=0.05,thresh=0.1,se.type="se.ols") {
    astd <- adjpow(tbd,pow,alpha,type="td",thresh=thresh,se.type=se.type)
    anls<-adjpow(tbd,pow,alpha,type="nls",thresh=thresh,se.type=se.type)
    anls2<-adjpow(tbd,pow,alpha,type="nls2",thresh=thresh,se.type=se.type)
    astd2<-adjpow(tbd,pow,alpha,type="td2",thresh=thresh,se.type=se.type)
    cbind(param,anls,anls2,astd,astd2)
}

findparam <- function(fun1,fun2,cf1,cf2,dk,method="BFGS") {
    mcf <- fun1(cf1,dk)
    fn <- function(p) {        
        sum((fun2(p,dk)-mcf)^2)
    }
    optim(cf2,fn,method=method,control=list(maxit=1000,reltol=sqrt(.Machine$double.eps)/10))
}

fp.2step <- function(fun1,fun2,cf1,cf2,dval) {
    res <- vector("list",length(dval))
    for(i in 1:length(dval)) {
        s1 <- findparam(fun1,fun2,cf1,cf2,dval[i],method="Nelder-Mead")
        s2 <- findparam(fun1,fun2,cf1,s1$par,dval[i],method="BFGS")
        res[[i]] <- list(s1=s1,s2=s2)
    }
    res
}
