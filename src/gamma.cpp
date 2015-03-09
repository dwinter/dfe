#include <cmath>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include <Rcpp.h>


struct convolve_params{
    double w;
    double Ve;
    double alpha;
    double beta;
    int* k;
};

double integrand(double x, void * p){
    //double w, double Ve, double alpha, double beta){
    convolve_params * params = (convolve_params *)p;
    double w = params->w;
    double sigma = sqrt(params->Ve);
    double k_alpha = params->alpha * (*params->k);
    double beta = params->beta;
    // gsl paramerization  is a=shape, b=scale (i.e. inverse of our rate)
    return(gsl_ran_gaussian_pdf((w-x), sigma) * gsl_ran_gamma_pdf(x, k_alpha,1/beta));
}

//' Density function of gamma dfe
//' @param obs, numeric, observed fitnesses
//' @param shape, numeric, shape of gamma distribution
//' @param rate, numeric, rate of gamma distribution
//' @param Ve, numeric, experimental variance
//' @param Ut, numeric, mutation rate
//' @param log, boolean, return log-liklihood (default=TRUE)
//' @export
// [[Rcpp::export]]
double dma_gamma(std::vector<double> obs, double shape, double rate, double Ve,  double Ut, bool log = true){
    double nobs = obs.size();
    int k = 1;
    double kfac = 1;
    std::vector<gsl_function> f_ptrs(nobs);
    std::vector<double> res (nobs);
    std::vector<convolve_params> param_vec (nobs);
    double running_prob = exp(-Ut);
    for(size_t i = 0; i < nobs ; i++){
        // create a vector of integrands
        gsl_function F;
        convolve_params p = {obs[i], Ve, shape, rate, &k};
        param_vec[i] = p;
        F.function = &integrand;
        F.params = &param_vec[i];
        f_ptrs[i] =  F;
        //Start each result with the k=0 case
        res[i] = gsl_ran_gaussian_pdf(obs[i], sqrt(Ve)) * running_prob;
    }
    double result;
    double error;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(10000);
    while (running_prob < 0.999){
        kfac *= k;
        double mu_prob = (exp(-Ut) * pow(Ut,k)) /kfac;
        for(size_t i = 0; i < nobs; ++i){
            gsl_integration_qagiu(&f_ptrs[i], 0., 1e-7, 1e-7, 10000, ws, &result, &error);
            res[i] += result * mu_prob;
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

