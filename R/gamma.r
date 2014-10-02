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

##'@export
dma_gamma <- function(w, s, B, Vc, Ut, log=FALSE){
    res <- sum(vapply(w, .dma_gamma, s=s, Vc=Vc, B=B, Ut=Ut, log=TRUE, FUN.VALUE=0.0))
    if(log){
        return(res)
    }
    return(exp(res))
}



##'@export
fit_ma_gamma <- function(obs, fixed=list(), start=list(), verbose=FALSE){
    all_args <- c("s", "B", "Vc", "Ut")
    known_args <- c( names(fixed), names(start) )  
    if(!all(all_args %in% known_args)){
       msg <- paste("Must set fixed or starting value for following params\n",
                    known_args[!(all_args %in% known_args )])
       stop(msg)
    }
    lower_bound <- c(s= 0, B=0, Ve=0, Ut=0)
    Q <- function(s, B, Vc, Ut){
        if(verbose){
            params <- match.call()
            print (sapply(as.list(params)[2:5], round, 4))
        }
        -dma_gamma(obs, s, B, Vc, Ut,log=TRUE)
    }
      
    mle(Q, start=start, fixed=fixed,
        method="L-BFGS-B", 
        lower=rep(1e-6,length(start)))
}



.dma_gamma <- function(w, s,B,Vc,Ut, log=FALSE){
    p_mu <- mu_scan(Ut)
    n <- length(p_mu) -1
    res <-  sum(sapply(0:n,  NG_convolution, z=w, s=s, B=B, Vc=Vc, verbose=FALSE) * p_mu)
    if(log){
        return(log(res))
    }
    return(res)
}
  

##' @export
moments_gamma <-function(s,B){
    return(c(mean=s/B, var=s/(B^2)))
}

NG_convolution <- function(z, s, Beta, Vc, k, verbose=FALSE){
    if(k==0){#gamma distr undefined
        return( dnorm(z, 0, sqrt(Vc)) )
    }
    integrand <- function(x,y){
        return(dnorm(y-x, 0, sqrt(Vc)) * dgamma(x, s*k, scale=Beta))
    }
    res <- integrate(integrand, z, lower= 0, upper=Inf)
    if(verbose){
        return(res)
    }
    else{
        return(res$value)
    }
}


