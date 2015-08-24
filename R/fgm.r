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
    if(length(Ve==1) & nloci > 1){
        idm <- diag(length(start))
        diag(idm) <- Ve
        Ve <- idm
    }
#    starting_fitness <- FUN(start)
    k <- rpois(n, Ut)
    experimental_variance <- rnorm(n, 0,Ve)

    mutate <- function(k){
        if(k==0){
           return(start)
        }
        else{
            mutations <- MASS::mvrnorm(k, rep(0, nloci), Ve)
            moves <- if(k==1) sum(mutations) else colSums(mutations)
            return(start + moves)
        }
    }
    new_positions <- vapply(k, mutate, FUN.VALUE=numeric(nloci))
    ending_fitnesses <- apply(new_positions, 2, FUN, start=start) + experimental_variance
    #return( ending_fitnesses - starting_fitness)
    -ending_fitnesses
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
estimate_FGM_moments <- function(n=1e5, start, Vm, FUN=sq_dist, plot=FALSE){
    sims <- replicate(n, start + rnorm(length(start), 0, sqrt(Vm)))
    w <-  apply(sims, 2, FUN) - FUN(start)
    mean_sim <- mean(w)
    var_sim <- var(w)
    if(plot){
        plot(density(w))
    }

    return(list(mean        =  mean_sim ,
                variance    =  var_sim, 
                skew_median =  (3*(mean_sim - median(w)))/ sqrt(var_sim),
                p_benifical = mean(w <0)))
    
}

