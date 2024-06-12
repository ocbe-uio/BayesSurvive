#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List func_MCMC_graph_cpp(
  Rcpp::List sobj,
  Rcpp::List hyperpar,
  Rcpp::List ini,
  uint S,
  std::string method,
  bool MRF_2b
) {
  Rcpp::List out = Rcpp::List::create();
  return out;
}
