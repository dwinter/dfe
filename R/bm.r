##' Bateman-mukai estimators for MA data
##'@export
##'@param delta_m numeric, change in mean fitness
##'@param delta_v numeric, change in among-line variance
##'@return U numeric BM estimator of mutation rate
##'@return Ea numeric BM estimator of effect-size
##'@examples
##' BM(-0.01, 0.001)
##' x <- rma_normal(n=100, s=0.01, Vs=1e-9, Ve=1e-3, Ut=10)
##' BM(mean(x), var(x) - 1e-3)

BM <- function(delta_m, delta_v){
    return(list(U=(2*delta_m^2)/delta_v, Ea = delta_v/delta_m))
}

