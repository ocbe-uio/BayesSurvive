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
  double tau = Rcpp::as<double>(hyperpar["tau"]);
  double cb = Rcpp::as<double>(hyperpar["cb"]);

  arma::vec beta_ini;
  arma::uvec acceptlee;

  if (method == "Pooled" && MRF_G) {
    arma::mat x = Rcpp::as<arma::mat>(sobj["X"]);
    uint J = Rcpp::as<uint>(hyperpar["J"]);
    arma::mat ind_r = Rcpp::as<arma::mat>(hyperpar["ind.r"]);
    arma::mat ind_d = Rcpp::as<arma::mat>(hyperpar["ind.d"]);
    arma::mat ind_r_d = Rcpp::as<arma::mat>(hyperpar["ind.r_d"]);
    double be_prop_sd_scale = Rcpp::as<double>(hyperpar["be.prop.sd.scale"]);
    arma::vec be_ini = Rcpp::as<arma::vec>(ini["beta.ini"]);
    arma::vec ga_ini = Rcpp::as<arma::vec>(ini["gamma.ini"]);
    arma::vec h = Rcpp::as<arma::vec>(ini["h"]);

    Rcpp::List erg = updateRP_genomic_cpp(
      p, x, J, ind_r, ind_d, ind_r_d, be_ini, be_prop_sd_scale, ga_ini, h, tau,
      cb
    );

    beta_ini = Rcpp::as<arma::vec>(erg["be.ini"]);
    acceptlee = Rcpp::as<arma::uvec>(erg["acceptl"]);
  } else {
    beta_ini = arma::zeros<arma::vec>(S);
    acceptlee = arma::zeros<arma::uvec>(S);
    for (uint g = 0; g < S; ++g) { // loop through subgroups
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
