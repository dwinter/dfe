#shamelessly stolen from RcppGSL

raise_missing_gsl <- function(){
    stop("Unable to detect compliation flags for GSL, is GSL installed?", call.=FALSE)
}

gsl_flags <- function(){
#    if(.Platform$OS.type == "windows"){
#        LIB_GSL <- Sys.getenv("LIB_GSL")
#        if(LIB_GSL == ""){
#            raise_missing_gsl()
#        }
#        return(sprintf( "-L%s/lib -lgsl -lgslcblas", LIB_GSL ))
    return(RcppGSL::LdFlags())
    }
    tryCatch(gsl_libs   <- system( "gsl-config --libs"   , intern = TRUE ), 
             error = raise_missing_gsl)
    gsl_libs
}

cat(gsl_flags())

     
   
