.libPaths("/scratch/lustre/zemlys/lib64/R/library")
library(midasr)
source("10func.R")
source("15shapes.R")
source("hahr.R")

xino <- function(n)rnorm(n)
library(snow)
library(Rmpi)
library(midasr)
library(MASS)
library(iterators)

print(sessionInfo())
cl<-makeMPIcluster(mpi.universe.size()-1)

ignore <- clusterEvalQ(cl, {
    library(midasr)
    library(foreach)
    library(Rsolnp)
    library(alabama)
    library(MASS)
    source("10func.R")
    source("15shapes.R")
    source("hahr.R")
    xino <- function(n)rnorm(n)    
NULL
})

load("param/hhp.100.RData")

oneh<-function(l) {
    gen <- gencr(l$N,l$n,
                 l$h0$wfun,l$h0$wpar,
                 l$h1$wfun,l$h1$wpar,
                 l$dk,xino,xino,                
                 m=l$m)
    save(gen,file=paste("results/hhp0617/res",l$i,l$n,l$cno,l$dk,l$m,"1000.RData",sep="-"))
    list(gen=gen,l=l)               
}

##Run comparison between hAhr and hAhr.fixed
tm<-snow.time(rr1<-clusterApplyLB(cl,hhp.100,oneh))
save(rr1,tm,file="hhp0617.RData")
