#include <RcppArmadillo.h>
#include "misc.h"
#include "updateRP_genomic_cpp.h"
// [[Rcpp::depends(RcppArmadillo)]]

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
Rcpp::List calJpost_cpp(
  const Rcpp::List sobj,
  const Rcpp::List hyperpar,
  const Rcpp::List ini,
  const unsigned int S,
  const std::string method,
  const bool MRF_G,
  const bool MRF_2b
) {
  // hyperparameters
  unsigned int p = Rcpp::as<unsigned int>(sobj["p"]);
  double c0 = Rcpp::as<double>(hyperpar["c0"]);
  double pi_ga = Rcpp::as<double>(hyperpar["pi.ga"]);
  double tau = Rcpp::as<double>(hyperpar["tau"]);
  double cb = Rcpp::as<double>(hyperpar["cb"]);

  if ((method == "CoxBVSSL" || method == "Sub-struct") || (method == "Pooled" && !MRF_G)) {
    if (method == "Pooled" && !MRF_G) {
      double lambda = Rcpp::as<double>(hyperpar["lambda"]);
      double pi_G = Rcpp::as<double>(hyperpar["pi.G"]);
    } else {
      double lambda = 0.0;
      double pi_G = 0.0;
    }
    double a = Rcpp::as<double>(hyperpar["a"]);
    double b = Rcpp::as<double>(hyperpar["b"]);
    arma::mat G_ini = Rcpp::as<arma::mat>(ini["G.ini"]);
  }

  double loglike, logpriorBeta, logpriorH, logpriorGamma, logjpost, logpriorOmega, logpriorX;
  if (method == "Pooled" && MRF_G) {
    Rcpp::List n = sobj["n"];
    arma::mat x = Rcpp::as<arma::mat>(sobj["X"]);
    arma::uvec J = arma::conv_to<arma::uvec>::from(list_to_vector(hyperpar["J"]));
    arma::mat ind_r_d = Rcpp::as<arma::mat>(hyperpar["ind.r_d"]);
    arma::mat ind_d = Rcpp::as<arma::mat>(hyperpar["ind.d"]);
    arma::vec hPriorSh = Rcpp::as<arma::vec>(hyperpar["hPriorSh"]);
    arma::vec beta_ini = Rcpp::as<arma::vec>(ini["beta.ini"]);
    arma::vec gamma_ini = Rcpp::as<arma::vec>(ini["gamma.ini"]);
    arma::vec h = Rcpp::as<arma::vec>(ini["h"]);
    arma::vec cbtau(p);
    for (uint i = 0; i < gamma_ini.n_elem; ++i) {
      cbtau(i) = tau * (gamma_ini(i) == 1 ? cb : 1);
    }

    Rcpp::List erg = calJpost_helper_cpp(cbtau, x, beta_ini, h, hPriorSh, c0, ind_r_d, ind_d);
    loglike = erg["loglike1"];
    logpriorBeta = erg["logpriorBeta1"];
    logpriorH = erg["logpriorH1"];
    logpriorGamma = arma::sum(gamma_ini * log(pi_ga)) + arma::sum((1 - gamma_ini) * log(1 - pi_ga));
    logjpost = loglike + logpriorGamma + logpriorBeta + logpriorH;
  } else {
    for (uint g = 0; g < S; ++g) {
      // n <- sobj$n[[g]]
      // X <- sobj$X[[g]]
      // J <- hyperpar$J[[g]]
      // ind.r_d <- hyperpar$ind.r_d[[g]]
      // ind.d <- hyperpar$ind.d[[g]]
      // hPriorSh <- hyperpar$hPriorSh[[g]]

      // gamma.ini <- ini$gamma.ini[[g]]
      // beta.ini <- ini$beta.ini[[g]]
      // h <- ini$h[[g]]
      // cbtau <- tau * ifelse(gamma.ini == 1, cb, 1)

      // erg <- calJpost_helper_cpp(cbtau, X, beta.ini, h, hPriorSh, c0, ind.r_d, ind.d)
      // loglike[g] <- erg$loglike1
      // logpriorBeta[g] <- erg$logpriorBeta1
      // logpriorH[g] <- erg$logpriorH1
      if (MRF_G) {
        // logpriorGamma[g] <- sum(gamma.ini * log(pi.ga)) + sum((1 - gamma.ini) * log(1 - pi.ga))
        // logjpost[g] <- loglike[g] + logpriorGamma[g] + logpriorBeta[g] + logpriorH[g]
      } else { // CoxBVSSL / Sub-struct model
        // C.ini <- ini$C.ini[[g]]
        // V.ini <- ini$V.ini[[g]]
        // Sig.ini <- ini$Sig.ini[[g]]
        // omega.mat <- matrix(0, p, p)
        for (uint i = 0; i < p; ++i) {
          // omega.mat[i, ] <- dnorm(C.ini[i, ], mean = rep(0, p), sd = sqrt(V.ini[i, ]), log = TRUE)
        }
        // diag(omega.mat) <- dexp(diag(C.ini), rate = lambda / 2, log = TRUE)

        // logpriorOmega[g] <- sum(omega.mat[upper.tri(omega.mat, diag = TRUE)])
        // logpriorX[g] <- sum(dmvnorm(X, mean = rep(0, p), sigma = Sig.ini, log = TRUE))
      }
    }
  }
  if (!MRF_G) {
    // pii.mat <- matrix(0, p * S, p * S)

    // id.mat is "full" graph (with all possible edges set to 1)
    if (method == "CoxBVSSL") {
      // id.mat <- do.call("cbind", rep(list(do.call("rbind", rep(list(diag(1, p, p)), S))), S))
    } else {
      // id.mat <- matrix(0, p * S, p * S)
    }
    for (uint g = 0; g < S; ++g) {
      // id.mat[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <- matrix(1, p, p)
    }
    // pii.mat[id.mat == 1 & G.ini == 0] <- log(1 - pi.G)
    // pii.mat[id.mat == 1 & G.ini == 1] <- log(pi.G)

    // logpriorGraph <- sum(pii.mat[upper.tri(pii.mat, diag = FALSE)])

    // gamma.vec <- unlist(ini$gamma.ini)

    if (MRF_2b) {
      for (uint g = 0; g < S; ++g) { // b1 * G_ss
        // G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <- b[1] * G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)]
      }
      for (uint g = 0; g < S - 1; ++g) { // b2 * G_rs
        for (uint r = g; r < S - 1; ++r) {
          // G.ini[(g - 1) * p + (1:p), r * p + (1:p)] <- G.ini[r * p + (1:p), (g - 1) * p + (1:p)] <- b[2] * G.ini[r * p + (1:p), (g - 1) * p + (1:p)]
        }
      }
      // b <- 1;
    }
    // logpriorGamma <- (a * sum(gamma.vec) + b * t(gamma.vec) %*% G.ini %*% gamma.vec)
    // logjpost <- sum(loglike) + sum(logpriorBeta) + sum(logpriorH) + logpriorGamma + sum(logpriorOmega) + sum(logpriorX) + logpriorGraph
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("loglike") = loglike, Rcpp::Named("logjpost") = logjpost
  );
  return out;
}
