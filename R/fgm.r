#' Simulate fitness effects under fishers geometric model 
#'
#' This function simulates fitness effects under a model in which the fitness
#' distribution is simulated as coming from Fisher's geometic model
#'
#'@export
#'@param n numeric, number of lines to simulate
#'@param start numeric,  vector describing the starting position relative to an
#' optimum at teh origin (can take any number of dimensions)
#'@param Vm numeric, variance in size of mutations
#'@param Ve numeric, envrionmental variance 
#'@param Ut numeric, expected number of mutations over the length of the
#' experiment
#'@return numeric, a vector of fitnesses
#' @param FUN, function, fitness funciton
#'@export
#'@examples
#' set.seed(123)
#' w <- rma_FGM(100, start=c(0.04,-0.04), Ve=1e-4, Vm=0.003, Ut=2)
#' mean(w)


rma_FGM <- function(n, start, Vm, Ve, Ut, FUN=gaussian_fitness){
    nloci <- length(start)
    if(length(Vm == 1) & nloci > 1){
        Vm <- make_sigma_mat(Vm, nloci)
    }
    k <- rpois(n, Ut)
    experimental_variance <- rnorm(n, 0, sqrt(Ve))
    new_positions <- vapply(k, mutate, start=start, Vm=Vm, FUN.VALUE=numeric(nloci))
    ending_fitnesses <- apply(new_positions, 2, FUN, start=start) + experimental_variance
    #return( ending_fitnesses - starting_fitness)
    -ending_fitnesses
}

make_sigma_mat <- function(Vm, nloci){
    idm <- diag(nloci)
    diag(idm) <- Vm
    Vm <- idm
    Vm
}

mutate <- function(start, Vm,k){
    if(k==0){
        return(start)
    } 
    mutations <- MASS::mvrnorm(k, start, sqrt(Vm))
    if(k >1){
       mutations <-  colSums(mutations)
    }
    mutations
}

gaussian_fitness <- function(coords, start, optimum="origin"){
    if(optimum=="origin"){
        optimum <- rep(0, length(coords))
    }
    Z <- dist(rbind(optimum, start))[1]
    dZ <- Z - dist(rbind(optimum, coords))[1]
    exp( -(Z-dZ)^2 / 2 ) - exp(-Z^2/2)
}



sq_dist <- function(coords, optimum="origin"){
    if(optimum=="origin"){
        optimum <- rep(0, length(coords))
    }
    d <- dist(rbind(optimum, coords))[1]
    return(d**2)
}



#' Estimate moments of d.f.e. as defined by fishers geometic model
#' @param n numeric number of observations from which to estimate the model
#' @param start, numeric vector of starting values
#' @param Vm, variance of mutation-size
#' @param FUN, function, fitness funciton
#' @param plot, boolean, make density plot of dfe (default=FALSE)
#' @export
estimate_FGM_moments <- function(n=1e4, start, Vm, FUN=gaussian_fitness, plot=TRUE){
    nloci <- length(start)
    if(length(Vm == 1) & nloci > 1){
        Vm <- make_sigma_mat(Vm, nloci)
    }
    k <- rep(1,n )
    new_positions <- vapply(k, mutate, start=start, Vm=Vm, FUN.VALUE=numeric(nloci))
    w <- -apply(new_positions, 2, FUN, start=start) 
    #return( ending_fitnesses - starting_fitness)
    mean_sim <- mean(w)
    var_sim <- var(w)
    den <- density(w)
    if(plot){
        plot(den)
    }
    return(list(mean        =  mean_sim ,
                variance    =  var_sim, 
                skew_median =  (3*(mean_sim - median(w)))/ sqrt(var_sim),
                p_benifical = mean(w <0),
                density     = den
                )
           )
    
}

