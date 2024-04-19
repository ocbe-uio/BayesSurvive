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

    double be_prop_me_ini = 0.;
    double be_prop_sd_ini = 0.;
    double D1 = 0.;
    double D2 = 0.;
    double D1_prop = 0.;
    double D2_prop = 0.;
    double loglh_ini = 0.;
    double loglh_prop = 0.;
    double logprior_prop = 0.;
    double logprior_ini = 0.;
    double logprop_prop = 0.;
    double logprop_ini = 0.;
    double logR = 0.;

    double be_prop_me = 1.;
    double be_prop_sd = 1.;

    arma::vec be_prop;
    arma::vec exp_xbeta;
    arma::mat h_exp_xbeta_mat, h_exp_xbeta_prop_mat;
    arma::mat exp_h_exp_xbeta_mat, exp_h_exp_xbeta_prop_mat; //, exp_xbeta_mat;
    arma::vec first_sum, second_sum;
    arma::vec first_sum_prop, second_sum_prop;
    arma::vec x_exp_xbeta, xbeta_prop, exp_xbeta_prop, x_exp_xbeta_prop, x_sq_exp_xbeta, x_sq_exp_xbeta_prop;
    arma::vec D1_1st, D1_2nd, D1_1st_prop, D1_2nd_prop;
    arma::vec D2_1st, D2_2nd, D2_1st_prop, D2_2nd_prop;
    arma::mat D1_2nd_den, D1_2nd_num, D1_2nd_den_prop, D1_2nd_num_prop;
    arma::mat D2_2nd_num, D2_2nd_den, D2_2nd_den_prop, D2_2nd_num_prop;

#endif
