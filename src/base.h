
#include <map>

typedef std::map<int, std::vector<double> > ProbMap;


ProbMap cache_mutation_probs(Rcpp::IntegerVector, double p_neutral);
