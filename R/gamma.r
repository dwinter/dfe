#' Simulate fitness effects under of Gamma model 
#'
#' This function simulates fitness effects under a model in which the fitness
#' distribution of mutations takes a Gamma distribution
#'
#'@export
#'@importFrom stats dpois
#'@param n numeric, number of lines to simulate
#'@param shape numeric,  shape parameter for Gamma
#'@param rate numeric, Scale paramater for Gamma
#'@param Ve numeric, envrionmental variance 
#'@param Ut numeric, expected number of mutations over the length of the
#' experiment
#'@return numeric, a vector of fitnesses
#'@examples
#' set.seed(123)
#' w <- rma_gamma(20, shape=0.1, rate=2, Ve=1e-4, Ut=3)
#' dma_gamma(obs=w, shape=0.1,rate= 2, Ve=1e-4, Ut = 3)
rma_gamma <- function(n, shape,rate, Ve, Ut){
    k <- rpois(n, Ut)
    between_line <- rnorm(n,0,sqrt(Ve))
    between_line +  rgamma(n, k*shape, rate)
}




#' Simulate fitness effects with a known number of mutations and Gamma DFE
#'
#' This function simulates fitness effects under a model in which the fitness
#' distribution of mutations takes a Gamma distribution
#'
#'@export
#'@importFrom stats rnorm rgamma median var rpois rbinom runif
#'@param n numeric, number of lines to simulate
#'@param shape numeric,  shape parameter for Gamma
#'@param rate numeric, Scale paramater for Gamma
#'@param Ve numeric, envrionmental variance 
#'@param k integer, total number of mutations in each line
#'@param p_neutral, proportion of mutations with no fitness effect
#'@return w, numeric simulate fitness of each line
#'@examples
#' k <- stats::rpois(20, 9)
#' w<- rma_known_gamma(shape=1, rate=20, Ve=0.01, k=k, p_neutral=0.4)
#' mean(w)

rma_known_gamma <- function(shape, rate, Ve, k, p_neutral){
    f <- function(m) {
        n <- length(m)
        between_line <- rnorm(n,0,sqrt(Ve))
        between_line +  rgamma(n, m*shape, rate)
    }
    rma_known_base(k, p_neutral, f)
}



# TODO
#Will have to over-write usage section for these
#
#' Fit a MA-model with Gamma dfe
#' @importFrom optimx optimx
#' @param obs, numeric observed fitnesses
#' @param shape numeric,  shape parameter for Gamma
#' @param rate numeric, Scale paramater for Gamma
#' @param Ve numeric, envrionmental variance 
#' @param Ut numeric, expected number of mutations over the length of the
#' @param fixed a named list containing values for any of the above paramaters,
#' that should be fixed during the maximum likelihood estimation
#' @param start a named list contaning starting values for any of the above
#' paramaters
#' Note, all paramaters must be given either a fixed or a starting value
#'@export


fit_ma_gamma <-  make_dfe_fitting_fxn(dma_gamma, 
                                      c("shape", "rate", "Ve", "Ut"), 
                                      lower = list(shape=1e-8, rate=1e-8, Ve=1e-8, Ut=0))


#' @export
fit_gamma_known <- make_dfe_fitting_fxn(dma_gamma_known, 
                                        c("shape", "rate", "Ve", "p_neutral", "k"),
                                        lower=list(shape=1e-8, rate=1e-8, Ve=1e-8, p_neutral=0),
                                        upper = list(p_neutral=1.0))


#fit_ma_gamma <- function(obs, fixed=list(), start=list(), verbose=FALSE){
#    all_args <- c("shape", "rate", "Ve", "Ut")
#    known_args <- c( names(fixed), names(start) )  
#    if(!all(all_args %in% known_args)){
#       msg <- paste("Must set fixed or starting value for following params\n",
#                    all_args[!(all_args %in% known_args )])
#       stop(msg)
#    }
#    lower_bound <- c(shape= 0, rate=0, Ve=0, Ut=0)
#    Q <- function(shape, rate, Ve, Ut){
#        if(verbose){
#            params <- match.call()
#            print (sapply(as.list(params)[2:5], round, 4))
#        }
#        if(any( c(shape, rate, Ve, Ut) <= 0)){
#           return(999999)
#        }
#        -dma_gamma(obs, shape, rate, Ve, Ut,log=TRUE)
#    }
#      
#    mle(Q, start=start, fixed=fixed,
#        method="L-BFGS-B",
#        lower=rep(1e-6,length(start)))
#        
#}
#
#' M.O.M estimate of gamma model with known mutation rate
#' @param obs, numeric, vector of observed fitnesses
#' @param Ve, numeric, experimental variance
#' @param Ut, numeric, mutation rate
#' @export
mom_ma_gamma <-function(obs, Ve, Ut){
    first <- mean(obs)
    second <- var(obs)
    denom <- first**2 - Ut * second   + Ut *Ve
    -(c(shape= first**2 / denom, rate=first*Ut/denom))
}



  
#' Get mean and variance of a Gamma distribution given shape and rate
#' paramaters
#' @param shape numeric,  shape parameter for Gamma
#' @param rate numeric, Scale paramater for Gamma
#' @export
moments_gamma <-function(shape, rate){
    return(c(mean=shape/rate, var=shape/(rate^2)))
}



