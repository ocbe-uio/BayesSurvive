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

    extern double be_prop_me_ini;
    extern double be_prop_sd_ini;
    extern double D1;
    extern double D2;
    extern double D1_prop;
    extern double D2_prop;
    extern double loglh_ini;
    extern double loglh_prop;
    extern double logprior_prop;
    extern double logprior_ini;
    extern double logprop_prop;
    extern double logprop_ini;
    extern double logR;

    extern double be_prop_me;
    extern double be_prop_sd;

    extern arma::vec be_prop;
    extern arma::vec exp_xbeta;
    extern arma::mat h_exp_xbeta_mat;
    extern arma::mat h_exp_xbeta_prop_mat;
    extern arma::mat exp_h_exp_xbeta_mat;
    extern arma::mat exp_h_exp_xbeta_prop_mat;
    extern arma::vec first_sum;
    extern arma::vec second_sum;
    extern arma::vec first_sum_prop;
    extern arma::vec second_sum_prop;
    extern arma::vec x_exp_xbeta;
    extern arma::vec xbeta_prop;
    extern arma::vec exp_xbeta_prop;
    extern arma::vec x_exp_xbeta_prop;
    extern arma::vec x_sq_exp_xbeta;
    extern arma::vec x_sq_exp_xbeta_prop;
    extern arma::vec D1_1st;
    extern arma::vec D1_2nd;
    extern arma::vec D1_1st_prop;
    extern arma::vec D1_2nd_prop;
    extern arma::vec D2_1st;
    extern arma::vec D2_2nd;
    extern arma::vec D2_1st_prop;
    extern arma::vec D2_2nd_prop;
    extern arma::mat D1_2nd_den;
    extern arma::mat D1_2nd_num;
    extern arma::mat D1_2nd_den_prop;
    extern arma::mat D1_2nd_num_prop;
    extern arma::mat D2_2nd_num;
    extern arma::mat D2_2nd_den;
    extern arma::mat D2_2nd_den_prop;
    extern arma::mat D2_2nd_num_prop;

#endif
