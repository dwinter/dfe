# Var[x] ~ mean/shape^3
# E[1/X] = 1/mean + 1/shape
#'@export
dinverse_gaussian <- function(x, mean, shape){
    sqrt(shape / (2*pi*x**3)) * exp( (-shape*(x-mean)**2)/(2*mean**2*x) )
}

#'@export
rinverse_gaussian <- function(n, mean, shape ){
    u <- mean
    y <- rnorm(n, 0, 1)**2
#    left <- mean + (mean*mean*y)/(2*shape) 
#    right <-  (mean/(2*shape)) * sqrt(4*mu*lambda*y + mu*mu*y*y)
    left <-  mean + (u**2*y)/(2*shape) 
    right <-  (mean / (2*shape)) * sqrt(4*mean*shape*y + mean**2*y**2)
    x <- left-right
    z <- runif(n, 0, 1)
    ifelse(z <= mean/(mean+x), x, mean**2/x)
}
 

rma_IG <- function(n, mean, shape, Ve, Ut){
    k <- rpois(n,Ut)
    between_line <- rnorm(n, 0, sqrt(Ve))
    mutations <- rinverse_gaussian(n, k*mean, k**2*shape)
    #deal with NAs cause by zero mutation cases
    ifelse(k > 0, between_line + mutations, between_line)    
}


#'@export
fit_ma_IG <- function(obs, fixed=list(), start=list(), verbose=FALSE){
    all_args <- list("mean", "shape", "Ve","Ut")
    known_args <- c(names(fixed), names(start))
    if(!all(all_args %in% known_args)){
        stop("Specify all the args yo")
    }
    Q <- function(mean, shape, Ve, Ut){
        if(any( c(shape, mean, Ve, Ut) <= 0)){
           return(999999)
        }
        res <- -dma_IG(obs, mean, shape, Ve, Ut,log=TRUE)
        if(verbose){
            params <- match.call()
            print (sapply(c(as.list(params)[2:5],res), round, 4))
        }
        res
    }
      
    mle(Q, start=start, fixed=fixed,
        method="L-BFGS-B",
        lower=rep(1e-6,length(start)))
        
}


#NIG_convolution <- function(z, mean, shape, Ve, k, verbose=FALSE){
#    if(k==0){
#        return(dnorm(z, 0, sqrt(Ve)))
#    }
#    res <- integrate(
#      function(x,y) dnorm(y-x, 0, sqrt(Ve)) * dinverse_gaussian(x, mean=mean*k, shape=shape*k**2),
#      y=z, lower=0, upper=Inf, abs.tol=1e-7
#    )
#    if(verbose){
#        return(res)
#    }
#    res$value
#}

