#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List func_MCMC_graph_cpp(
  const Rcpp::List sobj,
  const Rcpp::List hyperpar,
  const Rcpp::List ini,
  const uint S,
  const std::string method,
  const bool MRF_2b
) {
  // Extracting data
  Rcpp::List sobj_n = sobj["n"];
  int n = Rcpp::as<int>(sobj_n[0]);
  int p = Rcpp::as<int>(sobj["p"]);
  Rcpp::List SSig = sobj["SSig"];

  // Assembling output
  Rcpp::List Sig = Rcpp::List::create(); // TEMP: placeholder
  arma::mat G = arma::zeros<arma::mat>(p, p); // TEMP: placeholder
  Rcpp::List V = Rcpp::List::create(); // TEMP: placeholder
  Rcpp::List C = Rcpp::List::create(); // TEMP: placeholder

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("Sig.ini") = Sig,
    Rcpp::Named("G.ini") = G,
    Rcpp::Named("V.ini") = V,
    Rcpp::Named("C.ini") = C
  );
  return out;
}
