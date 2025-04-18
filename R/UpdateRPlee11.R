#' @title Update coefficients of Bayesian Cox Models
#'
#' @description
#' This an internal function to update coefficients of the Bayesian Cox Lasso Model
#'
#' @name UpdateRPlee11
#'
#' @inheritParams BayesSurvive
#' @inheritParams func_MCMC
#' @inheritParams func_MCMC_graph
#'
#' @return A list object with component 'beta.ini' for the updated coefficients
#' and component 'acceptlee' for the MCMC acceptance rate
#'
#' @export
UpdateRPlee11 <- function(sobj, hyperpar, ini, S, method, MRF_G, cpp = FALSE) {
  if (cpp) {
    return(UpdateRPlee11_cpp(sobj, hyperpar, ini, S, method, MRF_G))
  }
  p <- sobj$p
  tau <- hyperpar$tau
  cb <- hyperpar$cb

  if (method == "Pooled" && MRF_G) {
    x <- sobj$X
    J <- hyperpar$J
    ind.r <- hyperpar$ind.r
    ind.d <- hyperpar$ind.d
    ind.r_d <- hyperpar$ind.r_d
    be.prop.sd.scale <- hyperpar$be.prop.sd.scale
    be.ini <- ini$beta.ini
    ga.ini <- ini$gamma.ini
    h <- ini$h

    # erg = UpdateRP.lee11.helper(n, p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
    erg <- updateRP_genomic_cpp(
      p, x, J, ind.r, ind.d, ind.r_d,
      be.ini, be.prop.sd.scale, ga.ini, h, tau, cb
    )

    beta.ini <- as.vector(erg$be.ini)
    acceptlee <- erg$acceptl
  } else {
    beta.ini <- acceptlee <- vector("list", S)
    for (g in 1:S) { # loop through subgroups

      x <- sobj$X[[g]]
      J <- hyperpar$J[[g]]
      ind.r <- hyperpar$ind.r[[g]]
      ind.d <- hyperpar$ind.d[[g]]
      ind.r_d <- hyperpar$ind.r_d[[g]]
      be.ini <- ini$beta.ini[[g]]
      be.prop.sd.scale <- hyperpar$be.prop.sd.scale[[g]]
      ga.ini <- ini$gamma.ini[[g]]
      h <- ini$h[[g]]

      # erg = UpdateRP.lee11.helper(n, p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
      erg <- updateRP_genomic_cpp(
        p, x, J, ind.r, ind.d, ind.r_d,
        be.ini, be.prop.sd.scale, ga.ini, h, tau, cb
      )

      beta.ini[[g]] <- as.vector(erg$be.ini)
      acceptlee[[g]] <- erg$acceptl
    }
  }
  return(list(beta.ini = beta.ini, acceptlee = acceptlee))
}
