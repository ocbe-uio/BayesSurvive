#include <RcppArmadillo.h>
#include "misc.h"
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
  // p <- sobj$p
  // tau <- hyperpar$tau
  // cb <- hyperpar$cb

  // if (method == "Pooled" && MRF_G) {
  //   x <- sobj$X
  //   J <- hyperpar$J
  //   ind.r <- hyperpar$ind.r
  //   ind.d <- hyperpar$ind.d
  //   ind.r_d <- hyperpar$ind.r_d
  //   be.prop.sd.scale <- hyperpar$be.prop.sd.scale
  //   be.ini <- ini$beta.ini
  //   ga.ini <- ini$gamma.ini
  //   h <- ini$h

  //   # erg = UpdateRP.lee11.helper(n, p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
  //   erg <- updateRP_genomic_cpp(
  //     p, x, J, ind.r, ind.d, ind.r_d,
  //     be.ini, be.prop.sd.scale, ga.ini, h, tau, cb
  //   )

  //   beta.ini <- as.vector(erg$be.ini)
  //   acceptlee <- erg$acceptl
  // } else {
  //   beta.ini <- acceptlee <- vector("list", S)
  //   for (g in 1:S) { # loop through subgroups

  //     x <- sobj$X[[g]]
  //     J <- hyperpar$J[[g]]
  //     ind.r <- hyperpar$ind.r[[g]]
  //     ind.d <- hyperpar$ind.d[[g]]
  //     ind.r_d <- hyperpar$ind.r_d[[g]]
  //     be.ini <- ini$beta.ini[[g]]
  //     be.prop.sd.scale <- hyperpar$be.prop.sd.scale[[g]]
  //     ga.ini <- ini$gamma.ini[[g]]
  //     h <- ini$h[[g]]

  //     # erg = UpdateRP.lee11.helper(n, p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
  //     erg <- updateRP_genomic_cpp(
  //       p, x, J, ind.r, ind.d, ind.r_d,
  //       be.ini, be.prop.sd.scale, ga.ini, h, tau, cb
  //     )

  //     beta.ini[[g]] <- as.vector(erg$be.ini)
  //     acceptlee[[g]] <- erg$acceptl
  //   }
  // }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("beta.ini") = NA_REAL, // TEMP
    Rcpp::Named("acceptlee") = NA_REAL // TEMP
  );
  return out;
}
