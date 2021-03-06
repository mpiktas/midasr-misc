library(midasr)
library("multicore")
library("iterators")
library("foreach")
library("doRNG")
library("doMC")
registerDoMC(8)
source("code.R")

g0 <- c(10,2,-10)
g1 <- c(10,-10,-10)

param <- expand.grid(n=c(75,125,200,500,1000),d=c(11,23,47),m=c(12,24),ar=c(0.5,0.9))

param <- subset(param,!(m==12 & d==47))

set.seed(1234)

if(size) {
    tbd <- simtb(N,param,nealmon,g0,innov.sd=innov.sd,simplify=TRUE)
}else {
    dval <- unique(param[,"d"])
    paropt <- fp.2step(f.theta.kz3,nealmon,g1,g0,dval)
    h1st <- sapply(paropt,with,s2$par)
    colnames(h1st) <- as.character(dval)    
    ptb <- simtb(N,param,f.theta.kz3,g1,nealmon,h1st,innov.sd=innov.sd,simplify=TRUE)
}

