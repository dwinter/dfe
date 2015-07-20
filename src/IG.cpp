#include <Rcpp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <Rcpp.h>


using namespace Rcpp;

//Not expported as teh R-version is vectorized and fast enough for almost every
// imagainable case. Just cleaner to calculate individual values in C++ for other 
//functions
double dIG(double x, double mean, double shape){
    double result = sqrt(shape / (2*M_PI*pow(x,3))) * exp( (-shape*pow(x-mean,2))/(2*pow(mean,2)*x) ); 
    return(result);
}


struct ig_convolve_params{
    double w;
    double Ve;
    double mean;
    double shape;
    int* k;
};


double ig_integrand(double x, void * p){
    //double w, double Ve, double alpha, double beta){
    ig_convolve_params * params = (ig_convolve_params *)p;
    double w = params->w;
    double sigma = sqrt(params->Ve);
    double k_mean = params->mean * (*params->k);
    double k_shape = params->shape* pow((*params->k),2);
    // gsl paramerization  is a=shape, b=scale (i.e. inverse of our rate)
    return(gsl_ran_gaussian_pdf((w-x), sigma) * dIG(x, k_mean, k_shape));
}


//' Density function for inverse-gaussian dfe
//' @param obs numeric, observed fitnesses
//' @param mean numeric, mean effect-size of mutations
//' @param shape numeric, shape of underlying DFE
//' @param Ve numeric variance of experimental system
//' @param Ut numeric mutation rate
//' @param log boolean return log liklihood (default=TRUE)
//' @export
// [[Rcpp::export]]
double dma_IG(NumericVector obs, double mean, double shape, double Ve, double Ut, bool log = true){ 
    gsl_set_error_handler_off ();
    int nobs = obs.size();
    double running_prob = exp(-Ut);
    int k = 1;
    double kfac = 1;
    std::vector<gsl_function> f_ptrs(nobs);
    std::vector<double> res (nobs, 0.0);
    std::vector<ig_convolve_params> param_vec (nobs);
    for(size_t i = 0; i< nobs; i++){
        gsl_function F;
        ig_convolve_params p = {obs[i], Ve, mean, shape, &k};
        param_vec[i] = p;
        F.function = &ig_integrand;
        F.params = &param_vec[i];
        f_ptrs[i] =  F;

        //zero mutation case is just Noromal(0, Ve)
        res[i] = exp(-(pow(obs[i],2)/(2*Ve))) /  (sqrt(2*M_PI) * sqrt(Ve)) *  running_prob;
    }
    double result;
    double error;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(10000);
    while (running_prob < 0.999){
        kfac *= k;
        double mu_prob = (exp(-Ut) * pow(Ut,k)) /kfac;
        for(size_t i = 0; i < nobs; ++i){
            int err = 0;
            err = gsl_integration_qagiu(&f_ptrs[i], 0., 1e-7, 1e-7, 10000, ws, &result, &error);
            res[i] += result * mu_prob;
            if(err){
                Rf_warning("GSL intergration returned error %d", err);
//                std::cout << err << std::endl;
            }
        }
        running_prob += mu_prob;
        k += 1;
    }
    double final_result = 0;
    for(size_t i = 0; i < nobs; ++i){
        final_result += std::log(res[i]);
    }
    if(log){
        return(final_result);
    }
    return(exp(final_result));
}

