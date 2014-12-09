mu_scan <-function(Ut, tolerance=1e-35, max_mu=1e4){
    res <- dpois(0:max_mu, Ut)
    return (res[1: max(which(res > tolerance))])
}
