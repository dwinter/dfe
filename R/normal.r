
.dma_normal <- function(obs, s, Vs, Ve, Ut, log=FALSE){
    p_mu <- mu_scan(Ut)
    n <- length(p_mu) -1 
    res <- sum(dnorm(obs, (0:n)*s, sqrt( ((0:n)*Vs) + Ve)) * p_mu)
    if(log){
        return(log(res))
    }
    return(res)
}


##' Calculate log-likilhood for a single observation given 
##'
##' This function calculate the liklihood of a given observed fitness under
##' a model in which the fitness effects of mutations are normall distributed
##'
##'@export
##'@param obs numeric, a single observed fitness from a an MA lne
##'@param s numeric,  mean fitness of mutations
##'@param Vs numeric, variance in fitness of mutations
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


rma_normal <- function(n, s, Vs, Ve, Ut){
    k <- rpois(n, Ut)
    vapply(k, function(x) rnorm(1, x*s, sqrt((x*Vs)+Ve)), FUN.VALUE=0)
}

##' Find the maximum liklihood estimate of paramaters in MA model
##'@export
##'@param obs observed fitness values
##"@param fixed named list of model paramater values to fix
##'@param starts named list of starting values of varying parameters
##'@param verbose logical, be verbose
##'@return mle fit object
##'@examples
##' set.seed(123)
##' w <- rma_normal(20, 0.1, 0.01, 0.01, 1)
##' fit_ma_normal(w, fixed=list(Ve=0.01), starts=list(s=0.01))


fit_ma_normal <- function(obs, fixed=NULL, starts=NULL, verbose=TRUE){
    require(stats4)
    #Start by dealing with set/known/fixed arguments
    all_args <- c("s", "Vs", "Ve", "Ut")
    known_args <- c( names(fixed), names(starts) )
    to_set <- all_args[!(all_args %in% known_args )]
    if(length(to_set) > 0){
        if(verbose){
            cat("Setting starting values for following variables at random:\n")
            cat(to_set, "\n")
        }
        random_starts <- ifelse(to_set == "Ut", rpois(1,5), abs(rnorm(1,0,0.1)))
        names(random_starts) <- to_set
        starts <- c(starts, random_starts)
    }
    lower_bound <- c(s=-Inf, Vs=1e-5, Vc=1e-5, Ut=1e-5)
    Q <- function(s, Vs, Ve, Ut){
       if(verbose){
            params <- match.call()
            to_print <- sapply(as.list(params)[2:5], round, 4)
       }
       if(any(c(Vs,Ve,Ut) < 0)){
            return(.Machine$double.xmax)
       }
       res <- -dma_normal_cpp(obs, s, Vs, Ve, Ut, log=TRUE)
       if(verbose){
           print( c(to_print, "-LL"=round(res,2)))
       }
       res
    }
    mle(Q, start=starts, fixed=fixed, 
           method="L-BFGS-B", 
           lower=lower_bound[names(starts)],
           gr= function(s, Vs, Ve, Ut) numDeriv::grad(Q, c(s=s, Vs=Vs, Ve=Ve, Ut=Ut))
  #         lower=rep(0, 4),
   #        upper=rep(Inf,4))
   )
}

#' @useDynLib dfe
#' @importFrom Rcpp sourceCpp
#' @export
Rcpp::cppFunction('
double dma_normal_cpp(std::vector<double> obs, double a, double Va, double Ve, double Ut, bool log){ 
    //starting values for prob and res are for special case of k=0
    int n = obs.size();
    std::vector<double> res (n, 0.0);
    double running_prob = exp(-Ut);
    for(size_t i = 0; i< n; i++){
        res.push_back( (exp(-(pow(obs[i],2)/(2*Ve))) /  (sqrt(2*M_PI) * sqrt(Ve))) *  running_prob ) ;
    }
    uint64_t kfac = 1;
    uint16_t k= 1;
    while(running_prob < 0.9999){      
        kfac *= k;
        double mu_prob = (exp(-Ut) * pow(Ut,k)) /kfac;
        for(size_t i = 0; i < n; i++){
            res[i] += mu_prob * (  exp( -( pow(obs[i]-k*a,2)/(2*(Ve+ k*Va)) )) / (sqrt(2*M_PI) * sqrt(k * Va + Ve)) ) ;
        }
        k += 1;
        running_prob += mu_prob;
    }
    double lik = 0;
    for(size_t i = 0; i < n; i++){
        lik += std::log(res[i]);
    }
    if(log){
        return(lik);
    }
    return(exp(lik));
}')

f <- function(a, Va, Ut){
  res <- -dma_normal_cpp(w, a=a, Va=Va, Ve=1e-4, Ut=Ut, log=TRUE)
  print(paste("F:", a, Va, Ut, res))
  res
}

lazy_gr <- function(theta){
    print(theta)
    g <- function(theta) f(theta[1], theta[2], theta[3])
    numDeriv::grad(g, theta)
}

.grad_a <- function(w, a, Va, Ut){
   max_k <- length( dfe:::mu_scan(Ut, tolerance=1e-4)  )
sum(per_k(w, a, 0:max_k, Ut, Va, Ve=1e-4))
}

grad_a <- function(w, a, Va, Ut) sum(vapply(w, .grad_a, a=a, Va=Va, Ut=Ut, FUN.VALUE=0.0))

A <- function(w,a,k) w - a*k

B <- function(w, a, k, U, Va, Ve){
    exp(-U - (A(w,a,k)**2/ (2*(k*Va+Ve) ))) 
}
per_k <- function(w, a, k, U, Va, Ve){
   top <- B(w,a,k,U,Va,Ve) * k * U**k * A(w,a,k)
   bottom <- sqrt(2*pi) * (k*Va+Ve)**(3/2) * factorial(k)
   top / bottom
}

gradient_k <- function(k, w, a, Va, Ut, Ve=1e-4){
   #lots of resusable vars...
   A <-  w - a*k
   B <-   exp(-Ut - (A**2/ (2*(k*Va+Ve) ))) 
   tvar <- k*Va +Ve
   kfac <- factorial(k)
   ##Only need this once per loop!
#   divisor <- -dma_normal_cpp(w, a, Va, Ve, Ut, log=FALSE)
   root_2_pi <- sqrt(2*pi)
   dA <- (B * k * Ut**k * A) / (root_2_pi * tvar**(3/2) * kfac)

   dVa_left <- (B * k*Ut**k) / (2*root_2_pi * tvar**(3/2) * kfac)
   dVa_right<- (B * k*Ut**k) * A**2 / (2*root_2_pi * tvar**(5/2) * kfac)
   dVa<- (-dVa_left  + dVa_right)

   dUt_left  <- (B * k * Ut**(-1+k))/(root_2_pi * sqrt(tvar) * kfac)
   dUt_right <- (B  * Ut**k)/        (root_2_pi * sqrt(tvar) * kfac)
   dUt <- (dUt_left -  dUt_right)
   res <- c(dA, dVa, dUt)
   res

}


.grad <- function(w, a, Va, Ut, Ve=1e-4){
   max_k <- length( dfe:::mu_scan(Ut, tolerance=1e-4))
   rowSums(sapply(1:max_k, function(x) gradient_k(x, w, a, Va, Ut, Ve)))
}

grad <- function(theta){
    res <- sapply(w, function(x) .grad(x, theta[1], theta[2], theta[3], Ve=1e-4))
    rowSums(res)
}


g <- function(w, Ve, U){
    exp(-U - (w**2/(2*Ve))) / (sqrt(2*pi) * sqrt(Ve))
}

###top/-exp(dma <- normal <- cpp(w, a, Va, Ve, Utog=TRUE)) 
###stats4::mle(f, start=list(a=a/2, Va=Va, Ut=4), gr=grad)   
