mu_scan <-function(Ut, tolerance=1e-35, max_mu=1e4){
    res <- dpois(0:max_mu, Ut)
    return (res[1: max(which(res > tolerance))])
}


#For the random MA sims with known mutation. 
#This function takes a vecor of integers to treat as mutations counts, 
# a proporiton of those variants that are expected to be neutral and a a function 
# to generate gitness based on some DFE. The function is expected to take a vector of 
# mutation counts (where each mutation is one draw from teh DFE):
rma_known_base <- function(k, pn, f){
    with_effect <- rbinom(length(k), size=k, p=1-pn)
    f(with_effect)
}


