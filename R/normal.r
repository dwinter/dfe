
.dma_normal <- function(obs, s, Vs, Ve, Ut, log=FALSE){
    p_mu <- mu_scan(Ut)
    n <- length(p_mu) -1 
    res <- sum(dnorm(obs, (0:n)*s, sqrt( ((0:n)*Vs) + Ve)) * p_mu)
    if(log){
        return(log(res))
    }
    return(res)
}


##' Calculate log-likilhood for a single observation given 
##'
##' This function calculate the liklihood of a given observed fitness under
##' a model in which the fitness effects of mutations are normall distributed
##'
##'@export
##'@param obs numeric, a single observed fitness from a an MA lne
##'@param s numeric,  mean fitness of mutations
##'@param Vs numeric, variance in fitness of mutations
##'@param Ve numeric, envrionmental variance 
##'@param Ut numeric, expected number of mutations over the length of the
##' experiment
##'@param log logical, return log-transformed likilihod?
##'@return numeric, likilood of observed data given paramaters
##'@examples
##' Changinh mutation rates:
##' dma_normal(0.1, 0.1, 0.01, 0.01, 1)
##' dma_normal(0.1, 0.1, 0.01, 0.01, 3)



dma_normal <- function(obs, s, Vs, Ve, Ut, log=FALSE){
    res <- sum(vapply(obs, .dma_normal, s=s, Vs=Vs, Ve=Ve, Ut=Ut, log=TRUE, FUN.VALUE=0.0))
    if(!log){
       return(exp(res))
    }
    res
}


##' Simulate fitness effects under of normal model 
##'
##' This function simulates fitness effects under a model in which the fitness
##' distribution of mutations is normally distributed
##'
##'@export
##'@param n numeric, number of lines to simulate
##'@param s numeric,  mean fitness of mutations
##'@param Vs numeric, variance in fitness of mutations
##'@param Ve numeric, envrionmental variance 
##'@param Ut numeric, expected number of mutations over the length of the
##' experiment
##'@return numeric, a vector of fitnesses
##'@examples
##' set.seed(123)
##' w <- rma_normal(20, 0.1, 0.01, 0.01, 1)
##' dma_normal(w, 0.1, 0.01,0.01, 1)
##' dma_normal(w, 0.1, 0.01,0.01, 3)


rma_normal <- function(n, s, Vs, Ve, Ut){
    k <- rpois(n, Ut)
    vapply(k, function(x) rnorm(1, x*s, sqrt((x*Vs)+Ve)), FUN.VALUE=0)
}

##' Find the maximum liklihood estimate of paramaters in MA model
##'@export
##'@param obs observed fitness values
##"@param fixed named list of model paramater values to fix
##'@param starts named list of starting values of varying parameters
##'@param verbose logical, be verbose
##'@return mle fit object
##'@examples
##' w <- rma_normal(20, 0.1, 0.01, 0.01, 1)
##' fit_ma_normal(w, fixed=list(Ve=0.01), starts=list(s=0.01))


fit_ma_normal <- function(obs, fixed=NULL, starts=NULL, verbose=TRUE){
    require(stats4)
    #Start by dealing with set/known/fixed arguments
    all_args <- c("s", "Vs", "Ve", "Ut")
    known_args <- c( names(fixed), names(starts) )
    to_set <- all_args[!(all_args %in% known_args )]
    if(length(to_set) > 0){
        if(verbose){
            cat("Setting starting values for following variables at random:\n")
            cat(to_set, "\n")
        }
        random_starts <- ifelse(to_set == "Ut", rpois(1,5), abs(rnorm(1,0,0.1)))
        names(random_starts) <- to_set
        starts <- c(starts, random_starts)
    }
    lower_bound <- c(s=-Inf, Vs=1e-5, Vc=1e-5, Ut=1e-5)
    Q <- function(s, Vs, Ve, Ut){
        if(verbose){
            params <- match.call()
            print (sapply(as.list(params)[2:5], round, 4))
       }
        if(any(c(Vs,Ve,Ut) < 0)){
            return(999999999)
        }
       -dma_normal(obs, s, Vs, Ve, Ut, log=TRUE)
    }
    mle(Q, start=starts, fixed=fixed, 
#           method="L-BFGS-B", 

#           lower=lower_bound[names(starts)])
#           lower=rep(0, 4),
#           upper=rep(Inf,4))
        )
}










