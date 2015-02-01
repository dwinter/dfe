#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double dma_normal_cpp(NumericVector obs, double a, double Va, double Ve, double Ut, bool log){ 
    //starting values for prob and res are for special case of k=0
    int n = obs.size();
    std::vector<double> res (n, 0.0);
    double running_prob = exp(-Ut);
    for(size_t i = 0; i< n; i++){
        res[0] = exp(-(pow(obs[i],2)/(2*Ve))) /  (sqrt(2*M_PI) * sqrt(Ve)) *  running_prob;
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
}



// [[Rcpp::export]]
NumericVector grad_normal_cpp(NumericVector obs, double a, double Va, double Ve, double Ut, bool log){ 
    int n = obs.size();
    std::vector<double> dA (n, 0.0);
    std::vector<double> dV (n, 0.0);
    std::vector<double> dU (n, 0.0);
    double running_prob = exp(-Ut);
    double total_var = Ve;
    double two_root_pi = sqrt(2*M_PI);
    double expected_fitness = 0;
    for(size_t i = 0; i< n; i++){
        //Both dA and dV = 0 when k = 0, so we don't need to add them.
        // when k=0 dU  deptends only on U and Ve so do the relatively simple calc
        // without other variables
        dU[i] = -exp( -Ut - (pow(obs[i],2)/2*Ve)) / (two_root_pi * sqrt(Ve));
    }
    uint64_t kfac = 1;
    uint16_t k= 1;
    while(running_prob < 0.999999){      
        kfac *= k;
        running_prob += (exp(-Ut) * pow(Ut,k)) /kfac;
        total_var += Va;
        expected_fitness += a;
        for(size_t i = 0; i < n; i++){
            double A = obs[i] - expected_fitness ;
            double B = exp(-Ut - ( pow(A,2) / (2* total_var) ) );
            dA[i] += (B * k * pow(Ut,k) * A) / (two_root_pi * pow(total_var,1.5 ) * kfac);
            dV[i] +=  ( (B * pow(k*Ut,k)) * pow(A,2) / (2*two_root_pi * pow(total_var, 2.5) * kfac) ) - ( (B * k * pow(Ut,k)) / (2*two_root_pi * pow(total_var, 1.5) * kfac) ); 
            dU[i] += ( (B * k * pow(Ut, k-1)) / (two_root_pi * sqrt(total_var) * kfac) ) - ( (B * pow(Ut, k)) / (two_root_pi * sqrt(total_var) * kfac) ); 
        std::cout <<  (B * k * pow(Ut,k) * A) / (two_root_pi * pow(total_var,1.5 ) * kfac) << std::endl;
        }
        k += 1;
    }
    NumericVector res (3, 0.0);

    double divisor = 0;
    for(size_t i = 0; i < n; i++){
        divisor = dma_normal_cpp(obs, a, Va, Ve, Ut, false);
        res[0] += dA[i]/divisor;
        res[1] += dV[i]/divisor;
        res[2] += dU[i]/divisor;
    }
    return(res);
}



