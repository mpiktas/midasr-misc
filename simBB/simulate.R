library(midasr)
source("code.R")

g0 <- c(10,2,-10)

param <- expand.grid(n=c(75,125,200,500,1000),d=c(11,23,47),m=c(12,24),ar=c(0.5,0.9))

param <- subset(param,!(m==12 & d==47))

tbd <- foreach(i=1:1,.combine="c") %do% {
    n <- param$n[i]
    dk <- param$d[i]
    m <- param$m[i]
    ar <- param$ar[i]
    list(sim.rowdata(2000,n,dk,m,ar,weight=nealmon,cf=g0))
}
