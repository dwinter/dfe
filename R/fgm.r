##' Simulate fitness effects under fishers geometric model 
##'
##' This function simulates fitness effects under a model in which the fitness
##' distribution is simulated as coming from Fisher's geometic model
##'
##'@export
##'@param n numeric, number of lines to simulate
##'@param start numeric,  vector describing the starting position relative to an
##' optimum at teh origin (can take any number of dimensions)
##'@param Vm numeric, variance in size of mutations
##'@param Ve numeric, envrionmental variance 
##'@param Ut numeric, expected number of mutations over the length of the
##' experiment
##'@return numeric, a vector of fitnesses
##'@examples
##' set.seed(123)
##' w <- 



rma_FGM <- function(n, start, Vm, Ve, Ut, FUN=sq_dist){
    starting_fitness <- FUN(start)
    k <- rpois(n, Ut)
    experimental_variance <- rnorm(n, 0,Ve)
    mutate <- function(k){
        if(k==0){
           return(start)
        }
        else{
            return(rowSums(replicate(k, mutate_FGM(start, Vm))))
        }
    }
    new_positions <- vapply(k, mutate, FUN.VALUE=numeric(length(start)))
    ending_fitnesses <- vapply(new_positions, FUN, 0.0) + experimental_variance
    return( ending_fitnesses - starting_fitness)
}

mutate_FGM <- function(coords, V){
    coords + rnorm(length(coords), 0, sqrt(V))
}

sq_dist <- function(coords, optimum="origin"){
    if(optimum=="origin"){
        optimum <- rep(0, length(coords))
    }
    d <- dist(rbind(optimum, coords))[1]
    return(d**2)
}


estimate_FGM_moments <- function(n=1e5, start, Vm, FUN=sq_dist){
    sims <- replicate(n, mutate_FGM(start, Vm))
    w <-  vapply(sims, FUN, 0.0) - FUN(start)
    mean_sim <- mean(w)
    var_sim <- var(w)
    return(list(mean        =  mean_sim ,
                variance    =  var_sim, 
                skew_median =  (3*(mean_sim - median(w)))/ sqrt(var_sim)))
    
}

