hAhrfix.test <- function(x,PHI=vcovHAC(x$unrestricted,sandwich=FALSE)) {
    prep <- midasr:::prep_hAh(x)
    
    unrestricted <- x$unrestricted

    nyx <- nrow(x$model)
    nkx <- ncol(x$model)-1
    II <- diag(nkx)-prep$XtX %*% prep$Delta.0
    A0 <- nyx * ginv(t(prep$P)) %*% II %*% PHI %*% t(II) %*% ginv(prep$P)

    STATISTIC <- t(prep$h.0)%*%A0%*%prep$h.0
    
    names(STATISTIC) <- "hAhr"
    METHOD <- "hAh restriction test (robust version)"
    PARAMETER <- prep$dk-length(coef(x))
    PVAL <- 1-pchisq(STATISTIC,PARAMETER)
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD), 
        class = "htest")
}
