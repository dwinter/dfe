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
##'@param a numeric,  shape parameter for Gamma
##'@param B numeric, Scale paramater for Gamma
##'@param Ve numeric, envrionmental variance 
##'@param Ut numeric, expected number of mutations over the length of the
##' experiment
##'@return numeric, a vector of fitnesses
##'@examples
##' set.seed(123)
##' w <- rma_gamma(20, 0.1, 2, 1e-4, 3)
##' dma_gamma(w, 0.1, 2, 1e-4, 3)
rma_gamma <- function(n, shape,rate, Ve, Ut){
    k <- rpois(n, Ut)
    between_line <- rnorm(n,0,sqrt(Ve))
    between_line +  rgamma(n, k*shape, rate)
}


##' Calculate probability density of a set of fitness measures under a
##' Normal-Gamma convolution model
##'
##'@param a numeric,  shape parameter for Gamma
##'@param B numeric, Scale paramater for Gamma
##'@param Ve numeric, envrionmental variance 
##'@param Ut numeric, expected number of mutations over the length of the
##'@export
##'@examples
##' dma_gamma( c(0.1,0.01), 0.05, 2, 1e-4, 2)
##dma_gamma <- function(w, a, B, Ve, Ut, log=FALSE){
##    res <- sum(vapply(w, .dma_gamma, a=a, Ve=Ve, B=B, Ut=Ut, log=TRUE, FUN.VALUE=0.0))
##    if(log){
##        return(res)
##    }
##    return(exp(res))
##}


##'@param a numeric,  shape parameter for Gamma
##'@param B numeric, Scale paramater for Gamma
##'@param Ve numeric, envrionmental variance 
##'@param Ut numeric, expected number of mutations over the length of the
##'@param fixed a named list containing values for any of the above paramaters,
##' that should be fixed during the maximum likelihood estimation
##'@param start a named list contaning starting values for any of the above
##' paramaters
##' Note, all paramaters must be given either a fixed or a starting value
##'@export

fit_ma_gamma <- function(obs, fixed=list(), start=list(), verbose=FALSE){
    require(stats4)
    all_args <- c("a", "B", "Ve", "Ut")
    known_args <- c( names(fixed), names(start) )  
    if(!all(all_args %in% known_args)){
       msg <- paste("Must set fixed or starting value for following params\n",
                    known_args[!(all_args %in% known_args )])
       stop(msg)
    }
    lower_bound <- c(a= 0, B=0, Ve=0, Ut=0)
    Q <- function(a, B, Ve, Ut){
        if(verbose){
            params <- match.call()
            print (sapply(as.list(params)[2:5], round, 4))
        }
        if(any( c(a, B, Ve, Ut) <= 0)){
           return(999999)
        }
        -dma_gamma(obs, a, B, Ve, Ut,log=TRUE)
    }
      
    mle(Q, start=start, fixed=fixed)
        #method="L-BFGS-B,
        # ower=rep(1e-6,length(start)))
        
}



##.dma_gamma <- function(w, a,B,Ve,Ut, log=FALSE){
##    p_mu <- mu_scan(Ut)
##    n <- length(p_mu) -1
##    res <-  sum(sapply(0:n,  NG_convolution, z=w, a=a, B=B, Ve=Ve, verbose=FALSE) * p_mu)
##    if(log){
##        return(log(res))
##    }
##    return(res)
##}
  
##' Get mean and variance of a Gamma distribution given shape and rate
##' paramaters
##' @export
moments_gamma <-function(a,B){
    return(c(mean=a/B, var=a/(B^2)))
}


NG_convolution <- function(z, a, Beta, Ve, k, verbose=FALSE){
    if(k==0){#gamma distr undefined
        return( dnorm(z, 0, sqrt(Ve)) )
    }
    integrand <- function(x,y){
        return(dnorm(y-x, 0, sqrt(Ve)) * dgamma(x, shape=k*a, rate=Beta))
    }
    res <- integrate(integrand, z, lower= .Machine$double.xmin, upper=Inf,                     abs.tol=1e-7)
    if(verbose){
        return(res)
    }
    else{
        return(res$value)
    }
}


