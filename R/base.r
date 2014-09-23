mu_scan <-function(Ut, tolerance=1e-25, max_mu=1e3){
    res <- dpois(0:max_mu, Ut)
    return (res[1: max(which(res > tolerance))])
}
