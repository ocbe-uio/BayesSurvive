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
  double b = Rcpp::as<double>(hyperpar["b"]);

  arma::vec beta_ini = Rcpp::as<arma::vec>(ini["beta.ini"]);
  arma::vec gamma_ini = Rcpp::as<arma::vec>(ini["gamma.ini"]);

  // if (method %in% c("Pooled") && MRF_G) {
    // G.ini <- hyperpar$G
  // }

  // if (!MRF_G) {
    // G.ini <- ini$G.ini
  // }

  // two different b in MRF prior for subgraphs G_ss and G_rs
  // if (MRF_2b && !MRF_G) {
    // for (g in 1:S) { # b1*G_ss
      // G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <-
        // b[1] * G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)]
    // }
    // for (g in 1:(S - 1)) { # b2*G_rs
      // for (r in g:(S - 1)) {
        // G.ini[(g - 1) * p + (1:p), r * p + (1:p)] <-
          // G.ini[r * p + (1:p), (g - 1) * p + (1:p)] <-
          // b[2] * G.ini[r * p + (1:p), (g - 1) * p + (1:p)]
      // }
    // }
  // } else {
    // if (!MRF_G) {
      // G.ini <- G.ini * b
    // }
  // }

  // if (method == "Pooled" && MRF_G) {
    arma::vec post_gamma = arma::zeros<arma::vec>(p);

    // for (j in 1:p) {
      // beta <- beta.ini[j]

      // ga.prop1 <- ga.prop0 <- gamma.ini # gamma with gamma_g,j=1 or 0
      // ga.prop1[j] <- 1
      // ga.prop0[j] <- 0
      // ga.prop1 <- unlist(ga.prop1)
      // ga.prop0 <- unlist(ga.prop0)

      // wa <- (a * sum(ga.prop1) + t(ga.prop1) %*% G.ini %*% ga.prop1) +
        // dnorm(beta, mean = 0, sd = tau * cb, log = TRUE)
      // wb <- (a * sum(ga.prop0) + t(ga.prop0) %*% G.ini %*% ga.prop0) +
        // dnorm(beta, mean = 0, sd = tau, log = TRUE)

      // w_max <- max(wa, wb)
      // pg <- exp(wa - w_max) / (exp(wa - w_max) + exp(wb - w_max))

      // gamma.ini[j] <- as.numeric(runif(1) < pg)
      // post.gamma[j] <- pg
    // }
  // } else {
    // post.gamma <- rep(list(rep(0, p)), S)

    // if (MRF_G) {
      // for (g in 1:S) { # loop through subgroups
        // for (j in 1:p) {
          // wa <- dnorm((beta.ini[[g]])[j], mean = 0, sd = cb * tau) * pi
          // wb <- dnorm((beta.ini[[g]])[j], mean = 0, sd = tau) * (1 - pi)
          // pgam <- wa / (wa + wb)
          // u <- runif(1)
          // gamma.ini[[g]][j] <- ifelse(u < pgam, 1, 0)
          // post.gamma[[g]][j] <- pgam
        // }
      // }
    // } else { # CoxBVS-SL or Sub-struct model

      // for (g in 1:S) { # loop through subgroups
        // for (j in 1:p) {
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
        // }
      // }
    // }
  // }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("gamma.ini") = gamma_ini,
    Rcpp::Named("post.gamma") = post_gamma
  );
  return out;
}
