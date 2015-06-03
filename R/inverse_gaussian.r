# Var[x] ~ mean/shape^3
# E[1/X] = 1/mean + 1/shape
#'@export
dinverse_gaussian <- function(x, mean, shape){
    sqrt(shape / (2*pi*x**3)) * exp( (-shape*(x-mean)**2)/(2*mean**2*x) )
}

#'@export
rinverse_gaussian <- function(n, mean, shape){
    mean <- u
    y <- rnorm(n, 0, 1)**2
    left <-  u + (u**2*y)/(2*shape) 
    right <-  (u / 2*shape) * sqrt(4*u*shape*y + u**2*y**2)
    x <- left-right
    z <- runif(n, 0, 1)
    ifelse(z <= u/(u+x), x, u**2/x)
}
    

    


