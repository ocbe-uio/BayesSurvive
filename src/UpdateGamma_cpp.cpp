#include <RcppArmadillo.h>
#include "misc.h"
// [[Rcpp::depends(RcppArmadillo)]]

double calc_wa_wb(
  const arma::vec& ga,
  const arma::mat& G_ini,
  const double beta,
  const double a,
  const double stdev
) {
  double product = arma::as_scalar(a * arma::sum(ga) + ga.t() * G_ini * ga);
  double jitter = R::dnorm(beta, 0, stdev, true);
  return product + jitter;
}

double calc_pg(
  const arma::vec& ga_prop1,
  const arma::vec& ga_prop0,
  const arma::mat& G_ini,
  const double beta,
  const double a,
  const double tau,
  const double cb
) {
  double wa = calc_wa_wb(ga_prop1, G_ini, beta, a, tau * cb);
  double wb = calc_wa_wb(ga_prop0, G_ini, beta, a, tau);
  double w_max = std::max(wa, wb);
  return std::exp(wa - w_max) / (std::exp(wa - w_max) + std::exp(wb - w_max));
}

// [[Rcpp::export]]
Rcpp::List UpdateGamma_cpp(
  const Rcpp::List sobj,
  const Rcpp::List hyperpar,
  const Rcpp::List ini,
  const unsigned int S,
  const std::string method,
  const bool MRF_G,
  const bool MRF_2b
  ){
  // Update latent variable selection indicators gamma with either independent
  // Bernoulli prior (standard approaches) or with MRF prior.

  unsigned int p = Rcpp::as<unsigned int>(sobj["p"]);
  double tau = Rcpp::as<double>(hyperpar["tau"]);
  double cb = Rcpp::as<double>(hyperpar["cb"]);
  double pi = Rcpp::as<double>(hyperpar["pi.ga"]);
  double a = Rcpp::as<double>(hyperpar["a"]);
  arma::rowvec b = arma::rowvec(p, arma::fill::value(Rcpp::as<double>(hyperpar["b"])));

  arma::mat beta_ini = arma::zeros<arma::mat>(p, S);
  arma::mat gamma_ini = arma::zeros<arma::mat>(p, S);
  if (!Rf_isNewList(ini["beta.ini"])) {
    beta_ini.col(0) = Rcpp::as<arma::vec>(ini["beta.ini"]);
    gamma_ini.col(0) = Rcpp::as<arma::vec>(ini["gamma.ini"]);
  } else {
    beta_ini = list_to_matrix(ini["beta.ini"]);
    gamma_ini = list_to_matrix(ini["gamma.ini"]);
  }

  arma::mat G_ini = arma::zeros<arma::mat>(p * S, p * S);
  if (method == "Pooled" && MRF_G) {
    // G_ini is not needed if method != "Pooled" and MRF_G
    G_ini = Rcpp::as<arma::mat>(hyperpar["G"]);
  } else if (!MRF_G) {
    G_ini = Rcpp::as<arma::mat>(ini["G.ini"]);
  }

  // two different b in MRF prior for subgraphs G_ss and G_rs
  if (MRF_2b && !MRF_G) {
    // TODO: test this case. Not default!
    for (unsigned int g = 0; g < S; g++) {
      arma::uvec g_seq = arma::regspace<arma::uvec>(g * p, g * p + p - 1); // equivalent to (g - 1) * p + (1:p)
      G_ini.submat(g_seq, g_seq) *= b(0);
    }
    for (unsigned int g = 0; g < S - 1; g++) {
      for (unsigned int r = g; r < S - 1; r++) {
        arma::uvec g_seq = arma::regspace<arma::uvec>(g * p, g * p + p - 1); // equivalent to (g - 1) * p + (1:p)
        arma::uvec r_seq = arma::regspace<arma::uvec>(r * p + p, r * p + 2 * p - 1); // equivalent to r * p + (1:p)
        G_ini.submat(g_seq, r_seq) = b(1) * G_ini.submat(r_seq, g_seq);
        G_ini.submat(r_seq, g_seq) *= b(1);
      }
    }
  } else if (!MRF_G) {
    G_ini *= b(0);
  }

  arma::mat post_gamma = arma::zeros<arma::mat>(p, S);
  arma::mat ga_prop1 = arma::zeros<arma::mat>(p, S);
  arma::mat ga_prop0 = arma::zeros<arma::mat>(p, S);
  if (method == "Pooled" && MRF_G) {
    for (unsigned int j = 0; j < p; j++) {
      double beta = beta_ini(j);

      ga_prop1.col(0) = gamma_ini;
      ga_prop0.col(0) = gamma_ini;
      ga_prop1(j) = 1;
      ga_prop0(j) = 0;

      double pg = calc_pg(ga_prop1, ga_prop0, G_ini, beta, a, tau, cb);

      gamma_ini(j) = R::runif(0., 1.) < pg;
      post_gamma(j) = pg;
    }
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("gamma.ini") = gamma_ini,
      Rcpp::Named("post.gamma") = post_gamma
    );
    return out;
  } else {

    if (MRF_G) {
      for (unsigned int g = 0; g < S; g++) { // loop through subgroups
        for (unsigned int j = 0; j < p; j++) {
          double wa = R::dnorm(beta_ini(j, g), 0.0, tau * cb, true) * pi;
          double wb = R::dnorm(beta_ini(j, g), 0.0, tau, true) * (1. - pi);
          double pgam = wa / (wa + wb);
          double u = R::runif(0., 1.);
          gamma_ini(j, g) = u < pgam;
          post_gamma(j, g) = pgam;
        }
      }
    } else { // CoxBVS-SL or Sub-struct model

      for (unsigned int g = 0; g < S; g++) {
        for (unsigned int j = 0; j < p; j++) {
          double beta = beta_ini(j, g);

          ga_prop1 = gamma_ini;
          ga_prop0 = gamma_ini;
          ga_prop1(j, g) = 1;
          ga_prop0(j, g) = 0;

          double pg = calc_pg(arma::vectorise(ga_prop1), arma::vectorise(ga_prop0), G_ini, beta, a, tau, cb);
          gamma_ini(j, g) = R::runif(0., 1.) < pg;
          post_gamma(j, g) = pg;
        }
      }

    }

    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("gamma.ini") = gamma_ini,
      Rcpp::Named("post.gamma") = post_gamma
    );
    return out;
  }
}
