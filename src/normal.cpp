#include <Rcpp.h>

using namespace Rcpp;

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
}'



double dma_normal_cpp(std::vector<double> obs, double a, double Va, double Ve, double Ut, bool log){ 
    int n = obs.size();
    std::vector<double> dA (n, 0.0);
    std::vector<double> dV (n, 0.0);
    std::vector<double> dU (n, 0.0);
    double running_prob = exp(-Ut);
    
    for(size_t i = 0; i< n; i++){
        //speciac case k = 0, solve in mathmatica
    }
    uint64_t kfac = 1;
    uint16_t k= 0;
    while(running_prob < 0.9999){      



