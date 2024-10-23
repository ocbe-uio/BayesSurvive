#ifndef misc
#define misc

#include <RcppArmadillo.h>

// Function declaration
arma::vec list_to_vector(Rcpp::List r_list);
arma::mat list_to_matrix(Rcpp::List r_list);
arma::cube list_to_cube(Rcpp::List r_list);

#endif
