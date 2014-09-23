##NOTE: FOR CONSITENCY THIS FUNCTION SHOULD BE PARAMTERISED WITH
## MEAN AND VARIANCE OF FITNESS DISTRIBUTION IN FUTURE
##
##' Simulate fitness effects under of Gamma model 
##'
##' This function simulates fitness effects under a model in which the fitness
##' distribution of mutations takes a Gamma distribution
##'
##'@export
##'@param n numeric, number of lines to simulate
##'@param s numeric,  shape parameter for Gamma
##'@param B numeric, Scale paramater for Gamma
##'@param Ve numeric, envrionmental variance 
##'@param Ut numeric, expected number of mutations over the length of the
##' experiment
##'@return numeric, a vector of fitnesses
##'@examples
##' set.seed(123)
##' w <- rma_normal(20, 0.1, 1 0.01, 1)

rma_gamma <- function(n, s,B, Ve, Ut){
    res <- numeric(n)
    k <- rpois(n, Ut)
    between_line <- rnorm(n,0,Ve)
    for(i in 1:n){
        if(k[i] == 0){
            res[i] <- between_line[i]
    }   
        else {
            res[i] <- between_line[i] + sum(rgamma(k[i], shape=s,rate=B))
        }
    }
    return(res)
}

moments_gamma <-function(s,B){
    return(c(mean=s*B, var=s*B^2))
}
