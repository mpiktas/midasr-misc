library(midasr)
source("bbcode.R")
source("code.R")

g0 <- c(10,2,-10)
ar <- 0.5
m <- 12
n <- 75
k <- 0
dk <- m*k+m-3

pp <- function(p,dk)cumsum(nealmon(p,dk))

set.seed(1234)
Rprof()
system.time(rb <- f.eilute(1,n,pp,g0,pp,3,ar,sigma_u=1,m,k,dk))
Rprof(NULL)

set.seed(1234)
Rprof()
system.time(rz <- sim.rowdata2(1,n,dk-1,m,ar,weight0=nealmon,cf0=g0,simplify=FALSE))

Rprof(NULL)

set.seed(1234)
#mm <- mm.sim(n,dk,m,ar,nealmon(g0,dk))

