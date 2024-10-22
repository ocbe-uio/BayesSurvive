#include <RcppArmadillo.h>
#include "misc.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List UpdateRPlee11_cpp(
  const Rcpp::List sobj,
  const Rcpp::List hyperpar,
  const Rcpp::List ini,
  const uint S,
  const std::string method,
  const bool MRF_G
){
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("beta.ini") = NA_REAL, // TEMP
    Rcpp::Named("acceptlee") = NA_REAL // TEMP
  );
  return out;
}
