#include <RcppArmadillo.h>
#include "misc.h"
#include "updateRP_genomic_cpp.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List UpdateRPlee11_cpp(
  const Rcpp::List sobj,
  const Rcpp::List hyperpar,
  const Rcpp::List ini,
  const unsigned int S,
  const std::string method,
  const bool MRF_G
){
  unsigned int p = Rcpp::as<unsigned int>(sobj["p"]);
  double tau = Rcpp::as<double>(hyperpar["tau"]);
  double cb = Rcpp::as<double>(hyperpar["cb"]);

  arma::mat beta_ini(p, S, arma::fill::zeros);
  arma::umat acceptlee(p, S, arma::fill::zeros);
  arma::mat be_ini(p, S);
  arma::mat ga_ini(p, S);
  arma::uvec J(S);
  arma::vec be_prop_sd_scale(S);

  arma::field<arma::mat> ind_r(S);
  arma::field<arma::mat> ind_d(S);
  arma::field<arma::mat> ind_r_d(S);

  Rcpp::List erg;

  if (method == "Pooled" && MRF_G) {
    arma::field<arma::vec> h(S);
    unsigned int n = Rcpp::as<unsigned int>(sobj["n"]);
    arma::cube x(n, p, S, arma::fill::zeros);
    x.slice(0) = Rcpp::as<arma::mat>(sobj["X"]);
    be_ini.col(0) = Rcpp::as<arma::vec>(ini["beta.ini"]);
    ga_ini.col(0) = Rcpp::as<arma::vec>(ini["gamma.ini"]);

    J(0) = Rcpp::as<unsigned int>(hyperpar["J"]);
    ind_r(0) = Rcpp::as<arma::mat>(hyperpar["ind.r"]);
    ind_d(0) = Rcpp::as<arma::mat>(hyperpar["ind.d"]);
    ind_r_d(0) = Rcpp::as<arma::mat>(hyperpar["ind.r_d"]);
    be_prop_sd_scale(0) = Rcpp::as<double>(hyperpar["be.prop.sd.scale"]);
    double be_prop_sd_scale_value = be_prop_sd_scale(0); // Extract the first element
    h(0) = Rcpp::as<arma::vec>(ini["h"]);

    erg = updateRP_genomic_cpp(
      p, x.slice(0), J(0), ind_r(0), ind_d(0), ind_r_d(0),
      be_ini.col(0), be_prop_sd_scale_value, ga_ini.col(0),
      h(0), tau, cb
    );

    beta_ini.col(0) = Rcpp::as<arma::vec>(erg["be.ini"]);
    acceptlee.col(0) = Rcpp::as<arma::uvec>(erg["acceptl"]);
  } else {
    Rcpp::List h_list = ini["h"];
    Rcpp::List ind_r_list = hyperpar["ind.r"];
    Rcpp::List ind_d_list = hyperpar["ind.d"];
    Rcpp::List ind_r_d_list = hyperpar["ind.r_d"];
    arma::field<arma::vec> h(S);
    arma::field<arma::mat> ind_r(S);
    arma::field<arma::mat> ind_d(S);
    arma::field<arma::mat> ind_r_d(S);
    arma::cube x = list_to_cube(sobj["X"]);

    J = arma::conv_to<arma::uvec>::from(list_to_vector(hyperpar["J"]));
    be_ini = list_to_matrix(ini["beta.ini"]);
    be_prop_sd_scale = list_to_vector(hyperpar["be.prop.sd.scale"]);
    ga_ini = list_to_matrix(ini["gamma.ini"]);
    for (unsigned int g = 0; g < S; ++g) { // loop through subgroups
      double be_prop_sd_scale_value = be_prop_sd_scale(g);
      h(g) = Rcpp::as<arma::vec>(h_list[g]);
      ind_r(g) = Rcpp::as<arma::mat>(ind_r_list[g]);
      ind_d(g) = Rcpp::as<arma::mat>(ind_d_list[g]);
      ind_r_d(g) = Rcpp::as<arma::mat>(ind_r_d_list[g]);

      erg = updateRP_genomic_cpp(
        p, x.slice(g), J(g), ind_r(g), ind_d(g), ind_r_d(g),
        be_ini.col(g), be_prop_sd_scale_value, ga_ini.col(g),
        h(g), tau, cb
      );

      beta_ini.col(g) = Rcpp::as<arma::vec>(erg["be.ini"]);
      acceptlee.col(g) = Rcpp::as<arma::uvec>(erg["acceptl"]);
    }
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("beta.ini") = beta_ini,
    Rcpp::Named("acceptlee") = acceptlee
  );
  return out;
}
