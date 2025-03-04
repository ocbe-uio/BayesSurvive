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
  const uint S,
  const std::string method,
  const bool MRF_G,
  const bool MRF_2b
  ){
  // Update latent variable selection indicators gamma with either independent
  // Bernoulli prior (standard approaches) or with MRF prior.

  uint p = Rcpp::as<uint>(sobj["p"]);
  double tau = Rcpp::as<double>(hyperpar["tau"]);
  double cb = Rcpp::as<double>(hyperpar["cb"]);
  double pi = Rcpp::as<double>(hyperpar["pi.ga"]);
  double a = Rcpp::as<double>(hyperpar["a"]);
  arma::rowvec b = arma::rowvec(p, arma::fill::value(Rcpp::as<double>(hyperpar["b"])));

  arma::mat beta_ini = list_to_matrix(ini["beta.ini"]);
  arma::mat gamma_ini = list_to_matrix(ini["gamma.ini"]);

  arma::mat G_ini = arma::zeros<arma::mat>(p, p);
  if (method == "Pooled" && MRF_G) {
    // G_ini is not needed if method != "Pooled" and MRF_G
    G_ini = Rcpp::as<arma::mat>(hyperpar["G"]);
  } else if (!MRF_G) {
    G_ini = Rcpp::as<arma::mat>(ini["G.ini"]);
  }

  // two different b in MRF prior for subgraphs G_ss and G_rs
  if (MRF_2b && !MRF_G) {
    // TODO: test this case. Not default!
    for (arma::uword g = 0; g < S; g++) {
      arma::uvec g_seq = arma::regspace<arma::uvec>(g * p, g * p + p - 1); // equivalent to (g - 1) * p + (1:p)
      G_ini.submat(g_seq, g_seq) *= b(0);
    }
    for (arma::uword g = 0; g < S - 1; g++) {
      for (arma::uword r = g; r < S - 1; r++) {
        arma::uvec g_seq = arma::regspace<arma::uvec>(g * p, g * p + p - 1); // equivalent to (g - 1) * p + (1:p)
        arma::uvec r_seq = arma::regspace<arma::uvec>(r * p + p, r * p + 2 * p - 1); // equivalent to r * p + (1:p)
        G_ini.submat(g_seq, r_seq) = b(1) * G_ini.submat(r_seq, g_seq);
        G_ini.submat(r_seq, g_seq) *= b(1);
      }
    }
  } else if (!MRF_G) {
    G_ini *= b(0);
  }

  arma::mat post_gamma(S, p, arma::fill::zeros);
  if (method == "Pooled" && MRF_G) {
    for (arma::uword j = 0; j < p; j++) {
      // FIXME: why are indices flipped here w.r.t. the other cases?
      double beta = beta_ini(0, j);

      arma::vec ga_prop1 = gamma_ini.t();
      arma::vec ga_prop0 = gamma_ini.t();
      ga_prop1(j) = 1;
      ga_prop0(j) = 0;

      double pg = calc_pg(ga_prop1, ga_prop0, G_ini, beta, a, tau, cb);

      gamma_ini(0, j) = R::runif(0, 1) < pg;
      post_gamma(0, j) = pg;
    }
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("gamma.ini") = gamma_ini,
      Rcpp::Named("post.gamma") = post_gamma
    );
    return out;
  } else {
    if (MRF_G) {
      for (arma::uword g = 0; g < S; g++) { // loop through subgroups
        for (arma::uword j = 0; j < p; j++) {
          double wa = R::dnorm(beta_ini(j, g), 0, tau * cb, true) * pi;
          double wb = R::dnorm(beta_ini(j, g), 0, tau, true) * (1 - pi);
          double pgam = wa / (wa + wb);
          double u = R::runif(0, 1);
          gamma_ini(j, g) = u < pgam;
          post_gamma(g, j) = pgam;
        }
      }
    } else { // CoxBVS-SL or Sub-struct model
      for (arma::uword g = 0; g < S; g++) {
        for (arma::uword j = 0; j < p; j++) {
          double beta = beta_ini(j, g);

          arma::vec ga_prop1 = gamma_ini;
          arma::vec ga_prop0 = gamma_ini;
          ga_prop1(j, g) = 1;
          ga_prop0(j, g) = 0;

          double pg = calc_pg(ga_prop1, ga_prop0, G_ini, beta, a, tau, cb);
          gamma_ini(j, g) = R::runif(0, 1) < pg;
          post_gamma(g, j) = pg;
        }
      }
    }
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("gamma.ini") = arma::trans(gamma_ini),
      Rcpp::Named("post.gamma") = arma::trans(post_gamma)
    );
    return out;
  }
}
