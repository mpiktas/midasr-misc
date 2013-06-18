hhp <- foreach(r=iter(expand.grid(1:4,c(300,1000,2000)),by="row"),.combine=c) %do% {
    i<-as.numeric(r[1,1])
    list(list(h0=list(wfun=shfun[[i]],wpar=shpar[[i]]),
              h1=list(wfun=shfname[[i]],wpar=shpar[[i]]),
              i=i,
              n=as.numeric(r[1,2]),
              dk=shdk[[i]],
              m=12
              ))
}

hhp.100 <- do.call("c",lapply(hhp,function(l,N=1000,sp=10)
                   lapply(1:sp,function(s){
                       l$N<- round(N/sp,0)
                       l$cno<-s
                       l})
                   ))

hhp5 <- foreach(r=iter(expand.grid(1:5,c(300,1000,2000)),by="row"),.combine=c) %do% {
    i<-as.numeric(r[1,1])
    list(list(h0=list(wfun=shfun[[i]],wpar=shpar[[i]]),
              h1=list(wfun=shfname[[i]],wpar=shpar[[i]]),
              i=i,
              n=as.numeric(r[1,2]),
              dk=shdk[[i]],
              m=12
              ))
}

hhp5.100 <- do.call("c",lapply(hhp5,function(l,N=1000,sp=10)
                   lapply(1:sp,function(s){
                       l$N<- round(N/sp,0)
                       l$cno<-s
                       l})
                   ))

h1start <- function(dk,wfun,wpar,wfun1,wpar1,...) {
    optf <- function(p) {
        sum((wfun(wpar,dk)-wfun1(p,dk))^2)
    }
    opt <- optim(wpar1,optf)
    oo <- optim(opt$par,optf,...)
    list(opt=oo,res=data.frame(x=1:dk,y=wfun(wpar,dk),y1=wfun1(oo$par,dk)))
}
