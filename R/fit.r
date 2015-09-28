#
#mom_and_ll <-  function(w, Ve, k, p_n){
#  mom_est <- mom_ma_gamma(w, Ve, mean ( k ) * (1-p_n)  )
#  shape_mom <- mom_est[["shape"]]
#  rate_mom <- mom_est[["rate"]]
#
#  LL <- dma_gamma_known(w, Ve=Ve, shape=shape_mom,  rate=rate_mom, p_neutral = p_n, k=k)
#  c( LL_mom=LL, p_n=p_n, shape_mom = shape_mom, rate_mom = rate_mom )
#}
#
#
#mom_and_maxamize <- function(w, Ve, k, p_n){
#    res <- mom_and_ll(w, Ve, k, p_n)
#    mod <- fit_gamma_known(w, 
#                           start=list(shape=res[["shape_mom"]], rate=res[["rate_mom"]]),
#                           fixed=list(Ve=Ve, p_neutral = p_n), k)
#    c(res, mod@coef, maxLL=-mod@min)
#                           
#}
#
#
#
#AIC <- function(ests)  6- 2 *ests[5,]
#ests <- sapply(Ut, mom_and_ll, w=w, Ve=0.001)



check_args <- function(required){
    parent_args <- sys.frame(sys.parent())
    all_args <- c(parent_args$start, parent_args$fixed)
    if( !all (required %in% names(all_args) )) {
        msg <- paste("Must set fixed or starting paramater values for",
                     all_args[!(all_args %in% required)])
        stop(msg, .call=FALSE)                      
    }
}


process_constraints <- function(starting_vals, constraint, lower=TRUE){
    bounds <- constraint[names(starting_vals)]
    NO_CONSTRAINT <- if(lower)  -Inf else Inf
    res <- vapply(bounds, function(x) if(is.null(x)) NO_CONSTRAINT else x, FUN.VALUE=0.0)
    if(all(is.infinite(res))){
        return(NULL)
    }
    unname(res)
}


# closures to generating fitting functions
make_dfe_fitting_fxn <- function(Q, all_params, verbose=TRUE, lower=list(), upper=list(), ...){
    f <- function(obs, start=list(), fixed=list(), verbose=versbose){
    
        check_args(all_params)
        #Set up the specific fitting function (with fixed values given and (to
        #match optimx) variables given as a single vector)
        Q_args <- formals(Q)
        Q_args$obs <- obs
        Q_args[names(fixed)] <- fixed
        Q_args[names(start)] <- NULL

        #Are we using contraints? If so, let's set them up 
        lower_bound <- process_constraints(start, lower, lower=TRUE)
        upper_bound <- process_constraints(start, upper, lower=FALSE)
                    
        LL <- function(x) -do.call(Q, c(x, Q_args))
        if(is.null(upper_bound)){
            if(is.null(lower_bound)){                
                return( optimx(unlist(start), LL, ...) )
            }
            return( optimx(par=unlist(start), fn=LL, lower = lower_bound, method = "L-BFGS-B", ...) )
        }
        if(is.null(lower_bound)){
            return( optimx(unlist(start), fn=LL, upper=upper_bound, "L-BFGS-B", ...) )
        }
        optimx(unlist(start), LL, lower=lower_bound, upper=upper_bound, "L-BFGS-B", ...)
    }
    f
}


#' @export
fit_gamma_known <- make_dfe_fitting_fxn(dma_gamma_known, 
                                        c("shape", "rate", "Ve", "p_neutral", "k"),
                                        lower=list(shape=1e-8, rate=1e-8, Ve=1e-8, p_neutral=0),
                                        upper = list(p_neutral=1.0))

    


#mu <- rpois(100,10)
#
#w <- rma_known_gamma(shape=1, rate=10, Ve=0.01, k=mu, p_neutral=0.6)
#mod <- fit_gamma_known(w, 
#                       fixed=list(Ve=0.01, shape=1, rate=10),
#                       start=list(p_neutral = 0.1),
#                       k=mu)
#
#pn <- seq(0.1, 0.9, 0.05)
#res <- vapply(pn, mom_and_ll, w=w, Ve=0.01, k=mu, FUN.VALUE=rep(0.01,4))
