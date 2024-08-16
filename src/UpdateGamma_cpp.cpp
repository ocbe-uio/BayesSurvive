#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
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
  arma::vec b = Rcpp::as<arma::vec>(hyperpar["b"]);

  arma::vec beta_ini(p);
  arma::vec gamma_ini(p);
  if (MRF_G) { // TODO: move this inside loops?
    arma::vec beta_ini = Rcpp::as<arma::vec>(ini["beta.ini"]);
    arma::vec gamma_ini = Rcpp::as<arma::vec>(ini["gamma.ini"]);
  } else {
    Rcpp::List beta_ini_list = Rcpp::as<Rcpp::List>(ini["beta.ini"]);
    beta_ini = Rcpp::as<arma::vec>(beta_ini_list[0]); // FIXME: should change with g, not fixed at 0
    Rcpp::List gamma_ini_list = Rcpp::as<Rcpp::List>(ini["gamma.ini"]);
    gamma_ini = Rcpp::as<arma::vec>(gamma_ini_list[0]); // FIXME: should change with g, not fixed at 0
  }

  arma::mat G_ini = arma::zeros<arma::mat>(p, p);
  if (method == "Pooled" && MRF_G) {
    G_ini = Rcpp::as<arma::mat>(hyperpar["G"]);
  } else if (!MRF_G) {
    G_ini = Rcpp::as<arma::mat>(ini["G.ini"]);
  }

  // two different b in MRF prior for subgraphs G_ss and G_rs
  if (MRF_2b && !MRF_G) {
    // TODO: test this case. Not default!
    for (arma::uword g = 0; g < S; g++) {
      arma::uvec g_seq = arma::regspace<arma::uvec>(g * p, g * p + p - 1); // equivalent to (g - 1) * p + (1:p)
      G_ini.submat(g_seq, g_seq) *= b[0];
    }
    for (arma::uword g = 0; g < S - 1; g++) {
      for (arma::uword r = g; r < S - 1; r++) {
        arma::uvec g_seq = arma::regspace<arma::uvec>(g * p, g * p + p - 1); // equivalent to (g - 1) * p + (1:p)
        arma::uvec r_seq = arma::regspace<arma::uvec>(r * p + p, r * p + 2 * p - 1); // equivalent to r * p + (1:p)
        G_ini.submat(g_seq, r_seq) = b[1] * G_ini.submat(r_seq, g_seq);
        G_ini.submat(r_seq, g_seq) *= b[1];
      }
    }
  } else if (!MRF_G) {
    G_ini = G_ini * b;
  }

  if (method == "Pooled" && MRF_G) {
    arma::vec post_gamma = arma::zeros<arma::vec>(p);
    for (arma::uword j = 0; j < p; j++) {
      double beta = beta_ini(j);

      arma::vec ga_prop1 = gamma_ini;
      arma::vec ga_prop0 = gamma_ini;
      ga_prop1(j) = 1;
      ga_prop0(j) = 0;

      // ga.prop1 <- unlist(ga.prop1)
      // ga.prop0 <- unlist(ga.prop0)

      double wa = arma::as_scalar((a * arma::sum(ga_prop1) + ga_prop1.t() * G_ini * ga_prop1) + R::dnorm(beta, 0, tau * cb, true));
      double wb = arma::as_scalar((a * arma::sum(ga_prop0) + ga_prop0.t() * G_ini * ga_prop0) + R::dnorm(beta, 0, tau, true));

      double w_max = std::max(wa, wb);
      double pg = std::exp(wa - w_max) / (std::exp(wa - w_max) + std::exp(wb - w_max));

      gamma_ini(j) = R::runif(0, 1) < pg;
      post_gamma(j) = pg;
    }
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("gamma.ini") = gamma_ini,
      Rcpp::Named("post.gamma") = post_gamma
    );
    return out;
  } else {
    // post.gamma <- rep(list(rep(0, p)), S)
    Rcpp::List post_gamma = Rcpp::List::create(S);

    if (MRF_G) {
      for (arma::uword g = 0; g < S; g++) { // loop through subgroups
        for (arma::uword j = 0; j < p; j++) {
          // wa <- dnorm((beta.ini[[g]])[j], mean = 0, sd = cb * tau) * pi
          // wb <- dnorm((beta.ini[[g]])[j], mean = 0, sd = tau) * (1 - pi)
          // pgam <- wa / (wa + wb)
          // u <- runif(1)
          // gamma.ini[[g]][j] <- ifelse(u < pgam, 1, 0)
          // post.gamma[[g]][j] <- pgam
        }
      }
    } else { // CoxBVS-SL or Sub-struct model
      for (arma::uword g = 0; g < S; g++) {
        for (arma::uword j = 0; j < p; j++) {
          // beta <- (beta.ini[[g]])[j]

          // ga.prop1 <- ga.prop0 <- gamma.ini # gamma with gamma_g,j=1 or 0
          // ga.prop1[[g]][j] <- 1
          // ga.prop0[[g]][j] <- 0
          // ga.prop1 <- unlist(ga.prop1)
          // ga.prop0 <- unlist(ga.prop0)

          // wa <- (a * sum(ga.prop1) + t(ga.prop1) %*% G.ini %*% ga.prop1) +
            // dnorm(beta, mean = 0, sd = tau * cb, log = TRUE)
          // wb <- (a * sum(ga.prop0) + t(ga.prop0) %*% G.ini %*% ga.prop0) +
            // dnorm(beta, mean = 0, sd = tau, log = TRUE)

          // w_max <- max(wa, wb)
          // pg <- exp(wa - w_max) / (exp(wa - w_max) + exp(wb - w_max))

          // gamma.ini[[g]][j] <- as.numeric(runif(1) < pg)
          // post.gamma[[g]][j] <- pg
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
