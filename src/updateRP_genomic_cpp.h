#ifndef PSBC_H
#define PSBC_H

#include <vector>
#include <sstream>
#include <random>

#ifdef CCODE
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

    arma::mat matProdVec(const arma::mat, const arma::vec);
    arma::vec sumMatProdVec(const arma::mat, const arma::vec);
    Rcpp::List updateRP_genomic_cpp(const unsigned int p,
                                    const arma::mat,
                                    const unsigned int,
                                    arma::mat,
                                    arma::mat,
                                    arma::mat,
                                    arma::vec,
                                    const double,
                                    arma::vec,
                                    arma::vec,
                                    const double,
                                    const double);

#endif
