#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "updateRP_genomic_cpp.h"

arma::mat matProdVec(const arma::mat x, const arma::vec y)
{
    // multiply (element-wise) a matrix to a expanded vector

    arma::mat mat_y = arma::zeros<arma::mat>(y.n_elem, x.n_cols);
    mat_y.each_col() = y;
    arma::mat spanMat = x % mat_y; // elementwise product
    return spanMat;
}

arma::vec sumMatProdVec(const arma::mat x, const arma::vec y)
{
    // compute "arma::sum( matProdVec( ind_r_d_, exp_xbeta ).t(), 1 );"

    arma::vec spanVec = arma::zeros(x.n_cols);
    for (unsigned int i = 0; i < x.n_cols; ++i)
        spanVec(i) = arma::dot(x.col(i), y);
    return spanVec;
}

// [[Rcpp::export]]
arma::vec updateBH_cpp(const arma::mat x_,
                       const arma::vec beta_,
                       const unsigned int J_,
                       const arma::mat ind_r_d_,
                       const arma::vec hPriorSh_,
                       const arma::vec d_,
                       const double c0_)
{
    // update cumulative baseline harzard
    // update the increment h_j in the cumulative baseline hazard in each interval

    arma::vec xbeta_ = x_ * beta_;
    xbeta_.elem(arma::find(xbeta_ > 700)).fill(700.);
    arma::vec h_rate = c0_ + sumMatProdVec(ind_r_d_, arma::exp( xbeta_ ));

    // arma::vec shape = hPriorSh_ + d_;
    arma::vec h_ = arma::zeros<arma::vec>(J_);
    for (unsigned int j = 0; j < J_; ++j)
        // double h_rate = c0_ + arma::dot( ind_r_d_.col(j), arma::exp( xbeta_ ) );
        // h_(j) = arma::randg(arma::distr_param(hPriorSh_(j) + d_(j), 1. / h_rate(j)));
        h_(j) = R::rgamma( hPriorSh_(j) + d_(j), 1. / h_rate(j) );

    return h_;
}

// [[Rcpp::export]]
Rcpp::List updateBH_list_cpp(const Rcpp::List x_,
                             const Rcpp::List beta_,
                             const Rcpp::List J_,
                             const Rcpp::List ind_r_d_,
                             const Rcpp::List hPriorSh_,
                             const Rcpp::List d_,
                             const double c0_)
{
    // update cumulative baseline harzard
    // update the increment h_j in the cumulative baseline hazard in each interval

    unsigned int S = J_.size();
    Rcpp::List h_(S);

    for (unsigned int g = 0; g < S; ++g)
    {
        arma::mat x_tmp = x_[g];
        arma::vec beta_tmp = beta_[g];
        unsigned int J_tmp = J_[g];
        arma::mat ind_r_d_tmp = ind_r_d_[g];
        arma::vec hPriorSh_tmp = hPriorSh_[g];
        arma::vec d_tmp = d_[g];

        arma::vec xbeta_ = x_tmp * beta_tmp;
        xbeta_.elem(arma::find(xbeta_ > 700)).fill(700.);
        arma::vec h_rate = c0_ + sumMatProdVec(ind_r_d_tmp, arma::exp( xbeta_ ));

        arma::vec h_tmp = arma::zeros<arma::vec>(J_tmp);
        for (unsigned int j = 0; j < J_tmp; ++j)
            h_tmp(j) = arma::randg(arma::distr_param(hPriorSh_tmp(j) + d_tmp(j), 1. / h_rate(j)));

        h_[g] = h_tmp;
    }

    return h_;
}

// [[Rcpp::export]]
Rcpp::List  calJpost_helper_cpp(const arma::vec cbtau,
                                const arma::mat x_,
                                const arma::vec beta_,
                                const arma::vec h_,
                                const arma::vec hPriorSh_,
                                const double c0_,
                                const arma::mat ind_r_d_,
                                const arma::mat ind_d_)
{
    // subfunction to update joint posterior distribution

    arma::vec xbeta_ = x_ * beta_;
    xbeta_.elem(arma::find(xbeta_ > 700)).fill(700.);
    arma::vec exp_xbeta = arma::exp(xbeta_);

    double first_sum_ini = arma::accu(-h_ % sumMatProdVec(ind_r_d_, exp_xbeta));

    arma::mat h_exp_xbeta_mat = -arma::kron(exp_xbeta, h_.t());
    h_exp_xbeta_mat.elem(arma::find(h_exp_xbeta_mat > -1.0e-7)).fill(-1.0e-7);
    h_exp_xbeta_mat = arma::log(1.0 - arma::exp(h_exp_xbeta_mat));
    // double second_sum_ini = arma::accu(arma::sum((h_exp_xbeta_mat % ind_d_).t(), 1));
    double second_sum_ini = arma::accu(h_exp_xbeta_mat % ind_d_);
    double loglike1 = first_sum_ini + second_sum_ini;

    double logpriorBeta1 = 0.;
    for (unsigned int j = 0; j < beta_.size(); ++j)
    {
        // logpriorBeta1 += arma::log_normpdf( beta_(j), 0.0, cbtau(j) );
        logpriorBeta1 += R::dnorm( beta_(j), 0.0, cbtau(j), true);
    }

    double logpriorH1 = 0.;
    for (unsigned int j = 0; j < h_.size(); ++j)
    {
        logpriorH1 += R::dgamma( h_(j), hPriorSh_(j), 1. / c0_, true );
    }

    return Rcpp::List::create(
                              Rcpp::Named("loglike1") = loglike1,
                              Rcpp::Named("logpriorBeta1") = logpriorBeta1,
                              Rcpp::Named("logpriorH1") = logpriorH1
                              );
}


// [[Rcpp::export]]
Rcpp::List updateRP_genomic_cpp(const unsigned int p,
                                const arma::mat x_,
                                const unsigned int J_,
                                arma::mat ind_r_,
                                arma::mat ind_d_,
                                arma::mat ind_r_d_,
                                arma::vec be_,
                                const double be_prop_sd_scale,
                                arma::vec ga_,
                                arma::vec h_,
                                const double tau,
                                const double cb)
{
    // update coefficients of genomic variables via a MH sampler

    arma::uvec updatej = arma::randperm(p);

    arma::vec xbeta_ = x_ * be_;
    arma::vec sd_be_ = arma::ones<arma::vec>(ga_.n_elem);
    sd_be_.elem(arma::find(ga_ == 1.)).fill(cb);
    sd_be_ = sd_be_ * tau;
    arma::uvec sampleRPg_accept_ = arma::zeros<arma::uvec>(p);

    unsigned int j = 0;
    for (unsigned int j_id = 0; j_id < p; ++j_id)
    {
        // declare local variables
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
        arma::mat h_exp_xbeta_mat;
        arma::mat h_exp_xbeta_prop_mat;
        arma::mat exp_h_exp_xbeta_mat;
        arma::mat exp_h_exp_xbeta_prop_mat;
        arma::vec first_sum;
        arma::vec second_sum;
        arma::vec first_sum_prop;
        arma::vec second_sum_prop;
        arma::vec x_exp_xbeta;
        arma::vec xbeta_prop;
        arma::vec exp_xbeta_prop;
        arma::vec x_exp_xbeta_prop;
        arma::vec x_sq_exp_xbeta;
        arma::vec x_sq_exp_xbeta_prop;
        arma::vec D1_1st;
        arma::vec D1_2nd;
        arma::vec D1_1st_prop;
        arma::vec D1_2nd_prop;
        arma::vec D2_1st;
        arma::vec D2_2nd;
        arma::vec D2_1st_prop;
        arma::vec D2_2nd_prop;
        arma::mat D1_2nd_den;
        arma::mat D1_2nd_num;
        arma::mat D1_2nd_den_prop;
        arma::mat D1_2nd_num_prop;
        arma::mat D2_2nd_num;
        arma::mat D2_2nd_den;
        arma::mat D2_2nd_den_prop;
        arma::mat D2_2nd_num_prop;

        j = updatej(j_id);
        xbeta_.elem(arma::find(xbeta_ > 700)).fill(700.);
        exp_xbeta = arma::exp(xbeta_);
        x_exp_xbeta = x_.col(j) % exp_xbeta;
        // D1_1st = - h_ % arma::sum( matProdVec( ind_r_d_, x_exp_xbeta ).t(), 1 );
        D1_1st = -h_ % sumMatProdVec(ind_r_d_, x_exp_xbeta);

        h_exp_xbeta_mat = -arma::kron(exp_xbeta, h_.t());
        h_exp_xbeta_mat.elem(arma::find(h_exp_xbeta_mat > -1.0e-7)).fill(-1.0e-7);
        exp_h_exp_xbeta_mat = arma::exp(h_exp_xbeta_mat);
        D1_2nd_den = 1. - exp_h_exp_xbeta_mat;
        D1_2nd_num = matProdVec(exp_h_exp_xbeta_mat, x_exp_xbeta);
        D1_2nd = h_ % arma::sum((D1_2nd_num / D1_2nd_den % ind_d_).t(), 1);
        D1 = arma::sum(D1_1st + D1_2nd) - 1. / sd_be_(j) / sd_be_(j) * be_(j);

        x_sq_exp_xbeta = x_.col(j) % x_.col(j) % exp_xbeta;
        // D2_1st = - h_ % arma::sum( matProdVec( ind_r_d_, x_sq_exp_xbeta ).t(), 1 );
        D2_1st = -h_ % sumMatProdVec(ind_r_d_, x_sq_exp_xbeta);
        D2_2nd_den = D1_2nd_den % D1_2nd_den;
        D2_2nd_num = matProdVec(exp_h_exp_xbeta_mat, x_sq_exp_xbeta) % (1. - exp_h_exp_xbeta_mat + h_exp_xbeta_mat);
        D2_2nd = h_ % arma::sum((D2_2nd_num / D2_2nd_den % ind_d_).t(), 1);
        D2 = arma::accu(D2_1st + D2_2nd) - 1. / sd_be_(j) / sd_be_(j);

        be_prop_me = be_(j) - D1 / D2;
        be_prop_sd = be_prop_sd_scale / sqrt(-D2);
        be_prop = be_;

        // genomic version:
        be_prop(j) = R::rnorm( be_prop_me, be_prop_sd );
        // be_prop(j) = arma::randn(arma::distr_param(be_prop_me, be_prop_sd));
        xbeta_prop = xbeta_ - x_.col(j) * be_(j) + x_.col(j) * be_prop(j);
        xbeta_prop.elem(arma::find(xbeta_prop > 700)).fill(700.);
        exp_xbeta_prop = arma::exp(xbeta_prop);
        x_exp_xbeta_prop = x_.col(j) % exp_xbeta_prop;
        // D1_1st_prop = - h_ % arma::sum( matProdVec( ind_r_d_, x_exp_xbeta_prop ).t(), 1 );
        D1_1st_prop = -h_ % sumMatProdVec(ind_r_d_, x_exp_xbeta_prop);

        h_exp_xbeta_prop_mat = -arma::kron(exp_xbeta_prop, h_.t());
        h_exp_xbeta_prop_mat.elem(arma::find(h_exp_xbeta_prop_mat > -1.0e-7)).fill(-1.0e-7);
        exp_h_exp_xbeta_prop_mat = arma::exp(h_exp_xbeta_prop_mat);
        D1_2nd_den_prop = 1. - exp_h_exp_xbeta_prop_mat;
        D1_2nd_num_prop = matProdVec(exp_h_exp_xbeta_prop_mat, x_exp_xbeta_prop);
        D1_2nd_prop = h_ % arma::sum((D1_2nd_num_prop / D1_2nd_den_prop % ind_d_).t(), 1);
        D1_prop = arma::accu(D1_1st_prop + D1_2nd_prop) - 1. / sd_be_(j) / sd_be_(j) * be_prop(j);

        x_sq_exp_xbeta_prop = x_.col(j) % x_.col(j) % exp_xbeta_prop;
        // D2_1st_prop = -h_ % arma::sum( matProdVec( ind_r_d_, x_sq_exp_xbeta_prop ).t(), 1);
        D2_1st_prop = -h_ % sumMatProdVec(ind_r_d_, x_sq_exp_xbeta_prop);
        D2_2nd_den_prop = D1_2nd_den_prop % D1_2nd_den_prop;
        D2_2nd_num_prop = matProdVec(exp_h_exp_xbeta_prop_mat, x_sq_exp_xbeta_prop) % (1. - exp_h_exp_xbeta_prop_mat + h_exp_xbeta_prop_mat);
        D2_2nd_prop = h_ % arma::sum((D2_2nd_num_prop / D2_2nd_den_prop, ind_d_).t(), 1);
        D2_prop = arma::accu(D2_1st_prop + D2_2nd_prop) - 1. / sd_be_(j) / sd_be_(j);
        be_prop_me_ini = be_prop(j) - D1_prop / D2_prop;
        be_prop_sd_ini = be_prop_sd_scale / sqrt(-D2_prop);

        // first_sum = arma::sum( matProdVec( ind_r_d_, exp_xbeta ).t(), 1 );
        first_sum = sumMatProdVec(ind_r_d_, exp_xbeta);
        second_sum = arma::sum((arma::log(D1_2nd_den) % ind_d_).t(), 1);

        loglh_ini = arma::accu(-h_ % first_sum + second_sum);
        // first_sum_prop = arma::sum( matProdVec( ind_r_d_, exp_xbeta_prop ).t(), 1) ;
        first_sum_prop = sumMatProdVec(ind_r_d_, exp_xbeta_prop);
        second_sum_prop = arma::sum((arma::log(D1_2nd_den_prop) % ind_d_).t(), 1);
        loglh_prop = arma::accu(-h_ % first_sum_prop + second_sum_prop);

        logprior_prop = R::dnorm( be_prop(j), 0.0, sd_be_(j), true);
        logprior_ini = R::dnorm( be_(j), 0.0, sd_be_(j), true);
        logprop_prop = R::dnorm( be_prop(j), be_prop_me_ini, be_prop_sd_ini, true);
        logprop_ini = R::dnorm( be_(j), be_prop_me, be_prop_sd, true);
        /*logprior_prop = arma::log_normpdf(be_prop(j), 0.0, sd_be_(j));
        logprior_ini = arma::log_normpdf(be_(j), 0.0, sd_be_(j));
        logprop_prop = arma::log_normpdf(be_prop(j), be_prop_me_ini, be_prop_sd_ini);
        logprop_ini = arma::log_normpdf(be_(j), be_prop_me, be_prop_sd);*/
        logR = loglh_prop - loglh_ini + logprior_prop - logprior_ini + logprop_ini - logprop_prop;

        if( log( R::runif(0., 1.) ) < logR )
        // if (log(arma::randu()) < logR)
        {
            be_(j) = be_prop(j);
            xbeta_ = xbeta_prop;
            sampleRPg_accept_(j) = sampleRPg_accept_(j) + 1;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("be.ini") = be_,
        Rcpp::Named("acceptl") = sampleRPg_accept_
      );
}
