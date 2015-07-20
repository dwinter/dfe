#include <cmath>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <Rcpp.h>
#include "base.h"

using namespace Rcpp;

ProbMap cache_mutation_probs(IntegerVector k, double p_neutral){
    ProbMap result;
    IntegerVector K = unique(k);
    for(size_t i = 0; i < K.size(); i ++){
        std::vector<double> probs(K[i]+1, 0.0);
        for(size_t j = 0; j < K[i]+1; j++){
            probs[j] = gsl_ran_binomial_pdf(j,  p_neutral, K[i]);
        }
        result[ K[i] ] = probs;
    }
    return result;
}
