##' Bateman-mukai estimators for MA data
##'@export
##'@param obs, a vector of ftiness values
##'@param Ve, an estimate of within-line ('environmental') variance
##'@return U numeric BM estimator of mutation rate
##'@return Ea numeric BM estimator of effect-size
##'@examples
##'
##' low_var <- rma_normal(n=100, s=0.01, Vs=1e-9, Ve=1e-4, Ut=10)
##' BM(low_var, 1e-4)
##' high_var <- rma_normal(n=100, s=0.01, Vs=1e-3, Ve=1e-4, Ut=10)
##' BM(high_var, 1e-4)
BM <- function(obs, Ve){
    delta_m <- mean(obs)
    delta_v <- var(obs) - Ve
    return(list(U=(delta_m^2)/delta_v, Ea = delta_v/delta_m))
}

