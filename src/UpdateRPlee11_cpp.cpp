#include <RcppArmadillo.h>
#include "updateRP_genomic_cpp.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List UpdateRPlee11_cpp(
  const Rcpp::List sobj,
  const Rcpp::List hyperpar,
  const Rcpp::List ini,
  const uint S,
  const std::string method,
  const bool MRF_G
){
  uint p = Rcpp::as<uint>(sobj["p"]);
  uint n = Rcpp::as<uint>(sobj["n"]);
  double tau = Rcpp::as<double>(hyperpar["tau"]);
  double cb = Rcpp::as<double>(hyperpar["cb"]);

  arma::mat beta_ini(p, S, arma::fill::zeros);
  arma::umat acceptlee(p, S, arma::fill::zeros);
  arma::cube x(n, p, S);
  arma::mat be_ini(p, S);
  arma::mat ga_ini(p, S);
  arma::uvec J(S);
  arma::vec be_prop_sd_scale(S);

  arma::field<arma::mat> ind_r(S);
  arma::field<arma::mat> ind_d(S);
  arma::field<arma::mat> ind_r_d(S);
  arma::field<arma::vec> h(S);

  if (method == "Pooled" && MRF_G) {
    x.slice(0) = Rcpp::as<arma::mat>(sobj["X"]);
    be_ini.col(0) = Rcpp::as<arma::vec>(ini["beta.ini"]);
    ga_ini.col(0) = Rcpp::as<arma::vec>(ini["gamma.ini"]);

    J(0) = Rcpp::as<uint>(hyperpar["J"]);
    ind_r(0) = Rcpp::as<arma::mat>(hyperpar["ind.r"]);
    ind_d(0) = Rcpp::as<arma::mat>(hyperpar["ind.d"]);
    ind_r_d(0) = Rcpp::as<arma::mat>(hyperpar["ind.r_d"]);
    be_prop_sd_scale(0) = Rcpp::as<double>(hyperpar["be.prop.sd.scale"]);
    double be_prop_sd_scale_value = be_prop_sd_scale(0); // Extract the first element
    h(0) = Rcpp::as<arma::vec>(ini["h"]);

    Rcpp::List erg = updateRP_genomic_cpp(
      p, x.slice(0), J(0), ind_r(0), ind_d(0), ind_r_d(0),
      be_ini.col(0), be_prop_sd_scale_value, ga_ini.col(0),
      h(0), tau, cb
    );

    beta_ini.col(0) = Rcpp::as<arma::vec>(erg["be.ini"]);
    acceptlee.col(0) = Rcpp::as<arma::uvec>(erg["acceptl"]);
  } else {
    for (uint g = 0; g < S; ++g) { // loop through subgroups
      // Basically the same as the above block, but with the g-th element of
      // each. Everything is a list :()

      // x <- sobj$X[[g]]
      // J <- hyperpar$J[[g]]
      // ind.r <- hyperpar$ind.r[[g]]
      // ind.d <- hyperpar$ind.d[[g]]
      // ind.r_d <- hyperpar$ind.r_d[[g]]
      // be.ini <- ini$beta.ini[[g]]
      // be.prop.sd.scale <- hyperpar$be.prop.sd.scale[[g]]
      // ga.ini <- ini$gamma.ini[[g]]
      // h <- ini$h[[g]]

      // erg <- updateRP_genomic_cpp(
      //   p, x, J, ind.r, ind.d, ind.r_d,
      //   be.ini, be.prop.sd.scale, ga.ini, h, tau, cb
      // )

      // beta.ini[[g]] <- as.vector(erg$be.ini)
      // acceptlee[[g]] <- erg$acceptl
    }
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("beta.ini") = beta_ini,
    Rcpp::Named("acceptlee") = acceptlee
  );
  return out;
}
