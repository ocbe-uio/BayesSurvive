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
  // #define ARMA_64BIT_WORD

  unsigned int p = Rcpp::as<unsigned int>(sobj["p"]);
  double tau = Rcpp::as<double>(hyperpar["tau"]);
  double cb = Rcpp::as<double>(hyperpar["cb"]);

  arma::mat beta_ini = arma::zeros<arma::mat>(p, S);
  arma::umat acceptlee = arma::zeros<arma::umat>(p, S);
  arma::mat be_ini(p, S);
  arma::mat ga_ini(p, S);
  arma::uvec J(S);
  arma::vec be_prop_sd_scale(S);

  // It is not optimal to define x as a cube, since n can vary with groups
  // TODO: remove the reallocation of other Rcpp::List into new variables here
  arma::cube x;
  unsigned int n;
  if (!Rf_isNewList(sobj["X"])) {
    n = Rcpp::as<unsigned int>(sobj["n"]);
    x = arma::zeros<arma::cube>(n, p, S);
    x.slice(0) = Rcpp::as<arma::mat>(sobj["X"]);
  } else {
    n = Rcpp::as<Rcpp::List>(sobj["n"])[0];
    x = arma::zeros<arma::cube>(n, p, S);
    x = list_to_cube(sobj["X"]);
  }

  Rcpp::List erg;

  if (method == "Pooled" && MRF_G) {
    be_ini.col(0) = Rcpp::as<arma::vec>(ini["beta.ini"]);
    ga_ini.col(0) = Rcpp::as<arma::vec>(ini["gamma.ini"]);

    J(0) = Rcpp::as<unsigned int>(hyperpar["J"]);
    be_prop_sd_scale(0) = Rcpp::as<double>(hyperpar["be.prop.sd.scale"]);
    double be_prop_sd_scale_value = be_prop_sd_scale(0); // Extract the first element

    erg = updateRP_genomic_cpp(
      p, x.slice(0), J(0),
      Rcpp::as<arma::mat>(hyperpar["ind.r"]),
      Rcpp::as<arma::mat>(hyperpar["ind.d"]),
      Rcpp::as<arma::mat>(hyperpar["ind.r_d"]),
      be_ini.col(0), be_prop_sd_scale_value, ga_ini.col(0),
      Rcpp::as<arma::vec>(ini["h"]),
      tau, cb
    );

    beta_ini.col(0) = Rcpp::as<arma::vec>(erg["be.ini"]);
    acceptlee.col(0) = Rcpp::as<arma::uvec>(erg["acceptl"]);
  } else {
    if (!Rf_isNewList(ini["beta.ini"])) {
      J(0) = Rcpp::as<unsigned int>(hyperpar["J"]);
      be_prop_sd_scale(0) = hyperpar["be.prop.sd.scale"];
      be_ini = Rcpp::as<arma::mat>(ini["beta.ini"]);
      ga_ini = Rcpp::as<arma::mat>(ini["gamma.ini"]);
    } else {
      J = arma::conv_to<arma::uvec>::from(list_to_vector(hyperpar["J"]));
      be_prop_sd_scale = list_to_vector(hyperpar["be.prop.sd.scale"]);
      be_ini = list_to_matrix(ini["beta.ini"]);
      ga_ini = list_to_matrix(ini["gamma.ini"]);
    }

    for (unsigned int g = 0; g < S; ++g) { // loop through subgroups
      double be_prop_sd_scale_value = be_prop_sd_scale(g);

      erg = updateRP_genomic_cpp(
        p, x.slice(g), J(g),
        Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(hyperpar["ind.r"])[g]),
        Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(hyperpar["ind.d"])[g]),
        Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(hyperpar["ind.r_d"])[g]),
        be_ini.col(g), be_prop_sd_scale_value, ga_ini.col(g),
        Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(ini["h"])[g]),
        tau, cb
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
