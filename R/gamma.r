

#' Simulate fitness effects under of Gamma model 
#'
#' This function simulates fitness effects under a model in which the fitness
#' distribution of mutations takes a Gamma distribution
#'
#'@export
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
#'@param n numeric, number of lines to simulate
#'@param shape numeric,  shape parameter for Gamma
#'@param rate numeric, Scale paramater for Gamma
#'@param Ve numeric, envrionmental variance 
#'@param k integer, total number of mutations in each line
#'@param p_neutral
#'@return w, numeric simulate fitness of each line
#'@examples
#' k <- rpois(20, 9)
#' w<- rma_known_gamma(shape=1, rate=20, Ve=0.01, k=k, p_neutral=0.4))
#' mean(w)

rma_known_gamma <- function(shape, rate, Ve, k, p_neutral){
    f <- function(m) {
        n <- length(m)
        between_line <- rnorm(n,0,sqrt(Ve))
        between_line +  rgamma(n, k*shape, rate)
    }
    rma_known_base(k, p_neutral, f)
}






#' Calculate probability density of a set of fitness measures under a
#' Normal-Gamma convolution model
#'@param w numeric, observed data
#'@param a numeric,  shape parameter for Gamma
#'@param B numeric, Scale paramater for Gamma
#'@param Ve numeric, envrionmental variance 
#'@param Ut numeric, expected number of mutations over the length of the
#'@param log boolean, return log liklihood
#'@export
#'@examples
#' dma_gamma( c(0.1,0.01), 0.05, 2, 1e-4, 2)
dma_gamma_r <- function(w, a, B, Ve, Ut, log=FALSE){
    res <- sum(vapply(w, .dma_gamma, a=a, Ve=Ve, B=B, Ut=Ut, log=TRUE, FUN.VALUE=0.0))
    if(log){
        return(res)
    }
    return(exp(res))
}


#' Fit a MA-model with Gamma dfe
#' @importFrom stats4 mle
#' @param obs, numeric observed fitnesses
#' @param verbose, boolean, print values at each execution (default TRUE)
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

fit_ma_gamma <- function(obs, fixed=list(), start=list(), verbose=FALSE){
    all_args <- c("shape", "rate", "Ve", "Ut")
    known_args <- c( names(fixed), names(start) )  
    if(!all(all_args %in% known_args)){
       msg <- paste("Must set fixed or starting value for following params\n",
                    all_args[!(all_args %in% known_args )])
       stop(msg)
    }
    lower_bound <- c(shape= 0, rate=0, Ve=0, Ut=0)
    Q <- function(shape, rate, Ve, Ut){
        if(verbose){
            params <- match.call()
            print (sapply(as.list(params)[2:5], round, 4))
        }
        if(any( c(shape, rate, Ve, Ut) <= 0)){
           return(999999)
        }
        -dma_gamma(obs, shape, rate, Ve, Ut,log=TRUE)
    }
      
    mle(Q, start=start, fixed=fixed,
        method="L-BFGS-B",
        lower=rep(1e-6,length(start)))
        
}

#' M.O.M estimate of gamma model with known mutation rate
#' @param obs, numeric, vector of observed fitnesses
#' @param Ve, numeric, experimental variance
#' @param Ut, numeric, mutation rate
#' @export
mom_ma_gamma <-function(obs, Ut, Ve){
    first <- mean(obs)
    second <- var(obs)
    denom <- first**2 - Ut * second   + Ut *Ve
    abs(c(shape= first**2 / denom, rate=first*Ut/denom))
}


.dma_gamma <- function(w, a,B,Ve,Ut, log=FALSE){
    p_mu <- mu_scan(Ut)
    n <- length(p_mu) -1
    res <-  sum(sapply(0:n,  NG_convolution, z=w, a=a, B=B, Ve=Ve, verbose=FALSE) * p_mu)
    if(log){
        return(log(res))
    }
    return(res)
}
  
#' Get mean and variance of a Gamma distribution given shape and rate
#' paramaters
#' @param shape numeric,  shape parameter for Gamma
#' @param rate numeric, Scale paramater for Gamma
#' @export
moments_gamma <-function(shape, rate){
    return(c(mean=shape/rate, var=shape/(rate^2)))
}


NG_convolution <- function(z, a, Beta, Ve, k, verbose=FALSE){
    if(k==0){#gamma distr undefined
        cat( paste(k, dnorm(z, 0, sqrt(Ve)), sep="\t" ), "\n")
        return( dnorm(z, 0, sqrt(Ve)) )
    }
    integrand <- function(x,y){
        return(dnorm(y-x, 0, sqrt(Ve)) * dgamma(x, shape=k*a, rate=Beta))
    }
    res <- integrate(integrand, z, lower= 0, upper=Inf,  abs.tol=1e-7)
    cat(paste(res$value, k, sep="\t"), "\n")
    if(verbose){
        return(res)
    }
    else{
        return(res$value)
    }
}


