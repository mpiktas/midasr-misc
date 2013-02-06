library(numDeriv)
library(midasr)
## gal kokie nors kiti paketai reikalingi

f.theta.al3 = function(gamma, dk) {
    i <- (1:dk)/100
    gamma[1] * exp(gamma[2] * i + gamma[3] * i^2)/sum(exp(gamma[2] * i + gamma[3] * 
        i^2))
}

f.theta.kz3 = function(p, dk) {
    i <- (1:dk)/100
    (p[1] * i) * exp(p[2] * i + p[3] * i^2)
}

# f.theta.al3.sum=function(gamma,dk) { cumsum(f.theta.al3(gamma, dk)) }

be.isskirciu = function(rinkinys) {
    # paprastai parametru nudivergavimo atvejais statistikos buna stipriai didesnes
    # uz statistikas, gautas sukonvergavimo atvejais. Taigi jei statistiku
    # rinkinyje rasim dideli tuscia tarpa, vadinasi uz to tarpo esancios
    # statistikos buvo gautos, kai parametrai nudivergavo.
    return(rinkinys)
    tarpas = 100
    rinkinys = sort(rinkinys, na.last = NA)  #isremoovina galimas NA reiksmes
    l = length(rinkinys)
    if(l>1) {
        gaps = rinkinys[2:l] - rinkinys[1:(l - 1)]
        if (all(gaps <= tarpas)) 
          return(rinkinys) else {
              last.valid = which(gaps > tarpas)[1]
          }
        rinkinys[1:last.valid]
    }
    else rinkinys
}

f.eilute = function(rep, n, f_dgp, gamma_dgp, f_used, r_f_used, ar.koef, sigma_u, 
    m, k,d_k) {
    library(midasr)
    library(numDeriv)
    rinkinys.rnls = numeric(rep)
    rinkinys.nls = numeric(rep)
    rinkinys.rnls_nls = numeric(rep)
    rinkinys.rnls_nls2 = numeric(rep)
    rinkinys.T_d = numeric(rep)
    for (i in 1:rep) {
        ## duomenys##
       # d_k = m * k + m
        n_d = n
        v = arima.sim(n_d * m + 10 * m, model = list(ar = ar.koef), sd = 1)
        x = cumsum(v)
        x = ts(x, freq = m)
        v = ts(v, freq = m)
        V = mls(v, 0:(d_k - 1), m)
        V = V[(length(V[, 1]) - n_d + 1):length(V[, 1]), ]
        VtV = t(V) %*% V
        X = mls(x, d_k, m)
        X = X[(length(X[, 1]) - n_d + 1):length(X[, 1]), ]
        Z = cbind(V, X)
        theta = f_dgp(gamma_dgp, d_k)
        theta_plius = c(theta, theta[length(theta)])
        y = Z %*% theta_plius + rnorm(n_d, sd = sigma_u)

        ## vertinimas## OLS##
        tetai = solve(t(Z) %*% Z) %*% t(Z) %*% y
        theta.ols = tetai[1:d_k]
        theta_d.ols = tetai[d_k + 1]
        lik = y - X * theta_d.ols
        s2 = sum((y - Z %*% tetai)^2)/(n_d - length(tetai))
        
        ## restricted nls##
        fn_ralt = function(gamma) {
            res = lik - V %*% f_used(gamma, d_k)
            sum(res^2) + 1e+07 * (f_used(gamma, d_k)[d_k] - theta_d.ols)^2
        }
        prad_alt = rnorm(r_f_used)
        opt_f = optim(prad_alt, fn_ralt, method = "BFGS", control = list(maxit = 1000, 
            reltol = sqrt(.Machine$double.eps)/10), hessian = T)
        gamma.rnls = opt_f$par
        theta.rnls = f_used(gamma.rnls, d_k)
        
        ## nls##
        fn0 = function(gamma) {
            res = lik - V %*% f_used(gamma, d_k)
            sum(res^2)
        }
        prad.nls = rnorm(r_f_used)
        opt_f <- optim(prad.nls, fn0, method = "BFGS", control = list(maxit = 1000, 
            reltol = sqrt(.Machine$double.eps)/10), hessian = T)
        gamma.nls = opt_f$par
        theta.nls = f_used(gamma.nls, d_k)

        ## plot## plot(theta_plius) ##tikri parametrai lines(tetai) ##OLS
        ## lines(theta.rnls, col=2) lines(theta.nls, col=3)
        
        ## statistika nrls##
        P = matrix(c(rep(0, d_k - 1), 1), nrow = 1)  #apribojimo matrica
        D.rnls = jacobian(f_used, gamma.rnls, dk = d_k)
        Delt.rnls = D.rnls %*% ginv(1/n_d * t(D.rnls) %*% VtV %*% D.rnls) %*% t(D.rnls)
        narys = Delt.rnls %*% t(P) %*% ginv(P %*% Delt.rnls %*% t(P)) %*% P %*% Delt.rnls
        Sigma.rnls = ginv(VtV/n_d) - Delt.rnls + narys
        h.rnls = sqrt(n_d) * (theta.ols - theta.rnls)/sqrt(s2)
        T.rnls = t(h.rnls) %*% ginv(Sigma.rnls) %*% h.rnls
        rinkinys.rnls[i] = T.rnls
        
        ## statistika nls##
        D.nls = jacobian(f_used, gamma.nls, dk = d_k)
        Delt.nls = D.nls %*% ginv(1/n_d * t(D.nls) %*% VtV %*% D.nls) %*% t(D.nls)
        Sigma.nls = ginv(VtV/n_d) - Delt.nls
        h.nls = sqrt(n_d) * (theta.ols - theta.nls)/sqrt(s2)
        T.nls = t(h.nls) %*% ginv(Sigma.nls) %*% h.nls
        rinkinys.nls[i] = T.nls
      #  browser()
        ## statistika T.rnls_nls
        h.rnls_nls = sqrt(n_d) * (theta.rnls - theta.nls)/sqrt(s2)
        Sigma.rnls_nls = Delt.rnls %*% t(P) %*% ginv(P %*% Delt.rnls %*% t(P)) %*% 
            P %*% Delt.rnls
        Sigma.rnls_nls2 = Delt.nls %*% t(P) %*% ginv(P %*% Delt.nls %*% t(P)) %*% 
            P %*% Delt.nls
        T.rnls_nls = t(h.rnls_nls) %*% ginv(Sigma.rnls_nls) %*% h.rnls_nls
        T.rnls_nls2 = t(h.rnls_nls) %*% ginv(Sigma.rnls_nls2) %*% h.rnls_nls
        rinkinys.rnls_nls[i] = T.rnls_nls
        rinkinys.rnls_nls2[i] = T.rnls_nls2
        
        ## T_d testas, (H0:f(gamma,d)=theta_d)## duom##
        V = mls(v, 0:(d_k - 2), m)
        V = V[(length(V[, 1]) - n_d + 1):length(V[, 1]), ]
        X = mls(x, d_k - 1, m)
        X = X[(length(X[, 1]) - n_d + 1):length(X[, 1]), ]
        Z = cbind(V, X)
        ## ols## ols theta_d##
        tetai = solve(t(Z) %*% Z) %*% t(Z) %*% y
        theta.ols = tetai[1:(d_k - 1)]
        theta_d.ols = tetai[d_k]
        lik = y - X * theta_d.ols
        s2 = sum((y - Z %*% tetai)^2)/(n_d - length(tetai))
        ## nls##
        fn0 = function(gamma) {
            res = lik - V %*% f_used(gamma, d_k)[1:(d_k - 1)]
            sum(res^2)
        }
        prad.nls = rnorm(r_f_used)
        opt_f <- optim(prad.nls, fn0, method = "BFGS", control = list(maxit = 1000, 
            reltol = sqrt(.Machine$double.eps)/10), hessian = T)
        gamma.nls = opt_f$par[1:r_f_used]
        theta.nls = f_used(gamma.nls, d_k)[1:(d_k - 1)]
        theta_d.nls = f_used(gamma.nls, d_k)[d_k]
        ## plot## plot(theta) ##tikri parametrai lines(tetai, type='l') ##OLS
        ## lines(c(theta.nls, theta_d.nls), col=2) statistika T_d##
        h_d = sqrt(n_d) * (theta_d.nls - theta_d.ols)/sqrt(s2)
        D_visa = jacobian(f_used, gamma.nls, dk = d_k)
        D_pirm = D_visa[1:(d_k - 1), ]
        d_pask = D_visa[d_k, ]
        Sigma_d = t(d_pask) %*% ginv(1/n_d * t(D_pirm) %*% t(V) %*% V %*% D_pirm) %*% 
            d_pask
        T_d = h_d^2/Sigma_d
        rinkinys.T_d[i] = T_d
        
    }  ##monte skliaustelio pabaiga

    ## density plot if(rep>1) { par(mfrow=c(2,2)) plot(density(rinkinys.rnls),
    ## main='rnls, df=d-r+1') lines(density(rchisq(100000, df=d_k-2)), col=2)
    ## plot(density(rinkinys.nls), main='nls, df=d-r') lines(density(rchisq(100000,
    ## df=d_k-3)), col=2) plot(density(rinkinys.rnls_nls), main='rnls-nls, df=1')
    ## lines(density(rchisq(100000, df=1)), col=2) plot(density(rinkinys.T_d),
    ## main='H0:f(,d)=theta_d, df=1') lines(density(rchisq(100000,df=1)), col=2) }
    
    
    ## NA rm## rinkinys.rnls=as.numeric(na.omit(rinkinys.rnls))
    ## l.rnls=length(rinkinys.rnls) rinkinys.nls=as.numeric(na.omit(rinkinys.nls))
    ## l.nls=length(rinkinys.nls)
    ## rinkinys.rnls_nls=as.numeric(na.omit(rinkinys.rnls_nls))
    ## l.rnls_nls=length(rinkinys.rnls_nls)
    ## rinkinys.rnls_nls2=as.numeric(na.omit(rinkinys.rnls_nls2))
    ## l.rnls_nls2=length(rinkinys.rnls_nls2)
    ## rinkinys.T_d=as.numeric(na.omit(rinkinys.T_d)) l.T_d=length(rinkinys.T_d)
    
    rinkinys.rnls = be.isskirciu(rinkinys.rnls)
    l.rnls = length(rinkinys.rnls)
    rinkinys.nls = be.isskirciu(rinkinys.nls)
    l.nls = length(rinkinys.nls)
    rinkinys.rnls_nls = be.isskirciu(rinkinys.rnls_nls)
    l.rnls_nls = length(rinkinys.rnls_nls)
    rinkinys.rnls_nls2 = be.isskirciu(rinkinys.rnls_nls2)
    l.rnls_nls2 = length(rinkinys.rnls_nls2)
    rinkinys.T_d = be.isskirciu(rinkinys.T_d)
    l.T_d = length(rinkinys.T_d)
    
    ## size/power
    galia.rnls = sum(rinkinys.rnls > qchisq(0.95, df = d_k - r_f_used + 1))/l.rnls
    galia.nls = sum(rinkinys.nls > qchisq(0.95, df = d_k - r_f_used))/l.nls
    galia.rnls_nls = sum(rinkinys.rnls_nls > qchisq(0.95, df = 1))/l.rnls_nls
    galia.rnls_nls2 = sum(rinkinys.rnls_nls2 > qchisq(0.95, df = 1))/l.rnls_nls2
    galia.T_d = sum(rinkinys.T_d > qchisq(0.95, df = 1))/l.T_d
    
    ## return##
    ats = c(galia.rnls, galia.nls, galia.rnls_nls, galia.rnls_nls2, galia.T_d)
    ats = t(as.matrix(ats))
    colnames(ats) = c("rnls", "nls", "rnls-nls(1)", "rnls-nls(2)", "H0:f(,d)=theta_d")
#    aa <- as.list(sys.frame(1))

    #return(list(ats,data=aa))
    ats
}


f.lentele = function(rep, f_dgp, gamma_f_dgp, f_used, r_f_used, sigma_u) {
    ar = c(0.5, 0.9)
    m = c(12, 24)
    k = c(0, 1)
    n = c(75, 125, 200, 500, 1000)
    l.ar = length(ar)
    l.m = length(m)
    l.k = length(k)
    l.n = length(n)
    n.row = l.ar * l.m * l.k * l.n
    n.col = 5
    
    lentele.pav = matrix(rep(0, n.row * 4), ncol = 4)
    colnames(lentele.pav) = c("ar", "m", "d", "n_d")
    lentele = matrix(rep(0, n.row * n.col), ncol = n.col)
    colnames(lentele) = c("rnls", "nls", "rnls-nls1", "rnls-nls2", "T_d")
    eilute = 0
    estimate.time = TRUE
    if (rep < 30) 
        estimate.time = FALSE
    if (estimate.time) {
        time.start = proc.time()
        for (ar.i in 1:l.ar) for (m.i in 1:l.m) for (k.i in 1:l.k) for (n.i in 1:l.n) {
            eilute = eilute + 1
            temp = f.eilute(2, n[n.i], f_dgp, gamma_f_dgp, f_used, r_f_used, ar[ar.i], 
                sigma_u, m[m.i], k[k.i])
        }
        time.end = proc.time()
        time.elapsed = (time.end - time.start)["elapsed"]
        time.estimate = time.elapsed * rep/2
        print(paste("estimated time left is", round(time.estimate/60, 2), "minutes"))
    }
    
    
    eilute = 0
    for (ar.i in 1:l.ar) for (m.i in 1:l.m) for (k.i in 1:l.k) for (n.i in 1:l.n) {
        eilute = eilute + 1
        print(paste("case:", "ar =", ar[ar.i], "m =", m[m.i], "d =", m[m.i] * k[k.i] + 
            m[m.i] - 1, "n =", n[n.i]))
        lentele.pav[eilute, ] = c(ar[ar.i], m[m.i], m[m.i] * k[k.i] + m[m.i] - 1, 
            n[n.i])
        temp = f.eilute(rep, n[n.i], f_dgp, gamma_f_dgp, f_used, r_f_used, ar[ar.i], 
            sigma_u, m[m.i], k[k.i])
        print(temp)
        lentele[eilute, ] = temp
    }
    
    return(cbind(lentele.pav, lentele))
} 
