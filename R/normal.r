
##' Calculate log-likilhood for a single observation given 
##'
##' This function calculate the liklihood of a given observed fitness under
##' a model in which the fitness effects of mutations are normall distributed
##'
##'@export
##'@param obs numeric, a single observed fitness from a an MA lne
##'@param a numeric,  mean fitness of mutations
##'@param Va numeric, variance in fitness of mutations
##'@param Ve numeric, envrionmental variance 
##'@param Ut numeric, expected number of mutations over the length of the
##' experiment
##'@param log logical, return log-transformed likilihod?
##'@return numeric, likilood of observed data given paramaters
##'@examples
##' Changinh mutation rates:
##' dma_normal(0.1, 0.1, 0.01, 0.01, 1)
##' dma_normal(0.1, 0.1, 0.01, 0.01, 3)



#dma_normal <- function(obs, s, Vs, Ve, Ut, log=FALSE){
#    res <- sum(vapply(obs, .dma_normal, s=s, Vs=Vs, Ve=Ve, Ut=Ut, log=TRUE, FUN.VALUE=0.0))
#    if(!log){
#       return(exp(res))
#    }
#    res
#}


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


rma_normal <- function(n, a, Va, Ve, Ut){
    k <- rpois(n, Ut)
    rnorm(n, k*a, sqrt(k*Va+Ve))
}



#' Method of moments estimators for the normally distributed DFE
#' @param obs observed fitness values
#' @param Ve estimate of btween-line (experimental) variance
#' @return a Mean mutational effect size
#' @return Va Variance of d.f.e
#' @return Ut expected number of mutaitons over the experiment
#' @export
#' @examples
#' set.seed(123)
#' w <- rma_normal(n=20, a=0.1, Va=0.01, Ve=0.01,Ut=1)
#' mom_ma_normal(w, 0.01)
mom_ma_normal <- function(obs, Ve){
    first  <- mean(obs)
    second <- var(obs) - Ve
    third  <- mean( (obs - mean(obs))^3)
    A <- sqrt(9 * second**2 - 8 * third * first)
    c(a  = ((3*second)  - A) / (4*first), 
      Va = (((second * A)/first) - (3*second**2/first) + (4*third) ) / (8*first),
      Ut = (first * A + 3 * second * first) / (2 * third) )
    
}



#' Find the maximum liklihood estimate of paramaters in MA model
#' @useDynLib dfe
#' @importFrom Rcpp, sourceCpp
#' @export
#' @param obs observed fitness values
#' @param fixed named list of model paramater values to fix
#' @param starts named list of starting values of varying parameters
#' @param verbose logical, be verbose
#' @return mle fit object
#' @examples
#' set.seed(123)
#' w <- rma_normal(20, 0.1, 0.01, 0.01, 1)
#' fit_ma_normal(w, fixed=list(Ve=0.01), starts=list(s=0.01))


fit_ma_normal <- function(obs, fixed=NULL, starts=NULL, verbose=TRUE){
    require(stats4)
    #Start by dealing with set/known/fixed arguments
    all_args <- c("a", "Va", "Ve", "Ut")
    known_args <- c( names(fixed), names(starts) )
    to_set <- all_args[!(all_args %in% known_args )]
    ##TODO BETTER STARTING VALUES
    if(length(to_set) > 0){
        if(verbose){
            cat("Setting starting values for following variables at random:\n")
            cat(to_set, "\n")
        }
        random_starts <- ifelse(to_set == "Ut", rpois(1,5), abs(rnorm(1,0,0.1)))
        names(random_starts) <- to_set
        starts <- c(starts, random_starts)
    }
    lower_bound <- c(a=-Inf, Va=1e-5, Ve=1e-5, Ut=1e-5)
    Q <- function(a, Va, Ve, Ut){
       if(verbose){
            params <- match.call()
            to_print <- sapply(as.list(params)[2:5], round, 4)
       }
       if(any(c(Va,Ve,Ut) < 0)){
            return(.Machine$double.xmax)
       }
       res <- -dma_normal(obs, a, Va, Ve, Ut, log=TRUE)
       if(verbose){
           print( c(to_print, "-LL"=round(res,2)))
       }
       res
    }
    mle(Q, start=starts, fixed=fixed, 
#           method="CG", 
#           lower=lower_bound[names(starts)],
           gr= function(par) grad_normal(obs=obs, a=par[1], Va=par[2], Ut=par[3], Ve=fixed$Ve)
#            gr = function(par) print(par)
  #         lower=rep(0, 4),
   #        upper=rep(Inf,4))
   )
}


