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
 
#'@export
rma_IG <- function(n, mean, shape, Ve, Ut){
    k <- rpois(n,Ut)
    between_line <- rnorm(n, 0, sqrt(Ve))
    mutations <- rinverse_gaussian(n, k*mean, k**2*shape)
    #deal with NAs cause by zero mutation cases
    ifelse(k > 0, between_line + mutations, between_line)    
}

#'@export 
rma_known_IG <- function(mean, shape, Ve, k, p_neutral){
    f <- function(m){
        err <- rnorm(length(m), 0, sqrt(Ve))
        mu <- rinverse_gaussian(length(m), mean*m, shape*m**2)
        ifelse(m > 0, err+mu, err)
    }
    rma_known_base(k, p_neutral, f)
}

#' @export 
mom_ma_IG <- function(w, Ve, Ut){
    first <- mean(w)
    second <- var(w) - Ve
    c(mean  = first/Ut, 
      shape = first^3 / (Ut * (Ut*second-first^2)) 
    )
}

#'@export
fit_ma_IG <- make_dfe_fitting_fxn(dma_IG, 
                                  c("shape", "mean", "Ve", "Ut"),
                                  lower=list(mean=1e-8, rate=1e-8, Ve=1e-8, Ut=0))


#'@export
fit_IG_known <- make_dfe_fitting_fxn(dma_IG_known, 
                                        c("shape", "mean", "Ve", "p_neutral", "k"),
                                        lower=list(mean=1e-8, rate=1e-8, Ve=1e-8, p_neutral=0),
                                        upper = list(p_neutral=1.0))


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
