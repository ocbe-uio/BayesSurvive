#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


Rcpp::List calJpost_cpp(
  const Rcpp::List sobj,
  const Rcpp::List hyperpar,
  const Rcpp::List ini,
  const unsigned int S,
  const std::string method,
  const bool MRF_G,
  const bool MRF_2b
) {
  //hyperparameters
  //p <- sobj$p
  //c0 <- hyperpar$c0
  //pi.ga <- hyperpar$pi.ga
  //tau <- hyperpar$tau
  //cb <- hyperpar$cb
  //
  //if (method %in% c("CoxBVSSL", "Sub-struct") ||
  //  (method == "Pooled" && !MRF_G)) {
  //  lambda <- hyperpar$lambda
  //  a <- hyperpar$a
  //  b <- hyperpar$b
  //  pi.G <- hyperpar$pi.G
  //  G.ini <- ini$G.ini
  //}
  //
  //if (method == "Pooled" && MRF_G) {
  //  n <- sobj$n
  //  X <- sobj$X
  //  J <- hyperpar$J
  //  ind.r_d <- hyperpar$ind.r_d
  //  ind.d <- hyperpar$ind.d
  //  hPriorSh <- hyperpar$hPriorSh
  //
  //  gamma.ini <- ini$gamma.ini
  //  beta.ini <- ini$beta.ini
  //  h <- ini$h
  //  cbtau <- tau * ifelse(gamma.ini == 1, cb, 1)
  //
  //  erg <- calJpost_helper_cpp(cbtau, X, beta.ini, h, hPriorSh, c0, ind.r_d, ind.d)
  //  loglike <- erg$loglike1
  //  logpriorBeta <- erg$logpriorBeta1
  //  logpriorH <- erg$logpriorH1
  //
  //  logpriorGamma <- sum(gamma.ini * log(pi.ga)) + sum((1 - gamma.ini) * log(1 - pi.ga))
  //  logjpost <- loglike + logpriorGamma + logpriorBeta + logpriorH
  //} else {
  //  loglike <- logpriorBeta <- logpriorH <- logpriorGamma <-
  //    logjpost <- logpriorOmega <- logpriorX <- numeric()
  //
  //  for (g in 1:S) {
  //    n <- sobj$n[[g]]
  //    X <- sobj$X[[g]]
  //    J <- hyperpar$J[[g]]
  //    ind.r_d <- hyperpar$ind.r_d[[g]]
  //    ind.d <- hyperpar$ind.d[[g]]
  //    hPriorSh <- hyperpar$hPriorSh[[g]]
  //
  //    gamma.ini <- ini$gamma.ini[[g]]
  //    beta.ini <- ini$beta.ini[[g]]
  //    h <- ini$h[[g]]
  //    cbtau <- tau * ifelse(gamma.ini == 1, cb, 1)
  //
  //    erg <- calJpost_helper_cpp(cbtau, X, beta.ini, h, hPriorSh, c0, ind.r_d, ind.d)
  //    loglike[g] <- erg$loglike1
  //    logpriorBeta[g] <- erg$logpriorBeta1
  //    logpriorH[g] <- erg$logpriorH1
  //
  //    if (MRF_G) {
  //      logpriorGamma[g] <- sum(gamma.ini * log(pi.ga)) +
  //        sum((1 - gamma.ini) * log(1 - pi.ga))
  //      logjpost[g] <- loglike[g] + logpriorGamma[g] + logpriorBeta[g] + logpriorH[g]
  //    } else { # CoxBVSSL/ Sub-struct model
  //
  //      C.ini <- ini$C.ini[[g]]
  //      V.ini <- ini$V.ini[[g]]
  //      Sig.ini <- ini$Sig.ini[[g]]
  //
  //      omega.mat <- matrix(0, p, p)
  //      for (i in 1:p) {
  //        omega.mat[i, ] <- dnorm(C.ini[i, ],
  //          mean = rep(0, p),
  //          sd = sqrt(V.ini[i, ]), log = TRUE
  //        )
  //      }
  //      diag(omega.mat) <- dexp(diag(C.ini), rate = lambda / 2, log = TRUE)
  //
  //      logpriorOmega[g] <- sum(omega.mat[upper.tri(omega.mat, diag = TRUE)])
  //      logpriorX[g] <- sum(dmvnorm(X,
  //        mean = rep(0, p),
  //        sigma = Sig.ini, log = TRUE
  //      ))
  //    }
  //  }
  //}
  //if (!MRF_G) {
  //  pii.mat <- matrix(0, p * S, p * S)
  //
  //  # id.mat is "full" graph (with all possible edges set to 1)
  //  if (method == "CoxBVSSL") {
  //    id.mat <- do.call("cbind", rep(list(
  //      do.call("rbind", rep(list(diag(1, p, p)), S))
  //    ), S))
  //  } else {
  //    id.mat <- matrix(0, p * S, p * S)
  //  }
  //  for (g in 1:S) {
  //    id.mat[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <- matrix(1, p, p)
  //  }
  //
  //  pii.mat[id.mat == 1 & G.ini == 0] <- log(1 - pi.G)
  //  pii.mat[id.mat == 1 & G.ini == 1] <- log(pi.G)
  //
  //  logpriorGraph <- sum(pii.mat[upper.tri(pii.mat, diag = FALSE)])
  //
  //  gamma.vec <- unlist(ini$gamma.ini)
  //
  //  if (MRF_2b) {
  //    for (g in 1:S) { # b1 * G_ss
  //      G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <-
  //        b[1] * G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)]
  //    }
  //    for (g in 1:(S - 1)) { # b2 * G_rs
  //      for (r in g:(S - 1)) {
  //        G.ini[(g - 1) * p + (1:p), r * p + (1:p)] <-
  //          G.ini[r * p + (1:p), (g - 1) * p + (1:p)] <-
  //          b[2] * G.ini[r * p + (1:p), (g - 1) * p + (1:p)]
  //      }
  //    }
  //    b <- 1
  //  }
  //  logpriorGamma <- (a * sum(gamma.vec) + b * t(gamma.vec) %*% G.ini %*% gamma.vec)
  //
  //  logjpost <- sum(loglike) + sum(logpriorBeta) + sum(logpriorH) +
  //    logpriorGamma + sum(logpriorOmega) + sum(logpriorX) + logpriorGraph
  //}
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("loglike") = Rcpp::NumericVector::create(),
    Rcpp::Named("logjpost") = Rcpp::NumericVector::create()
  );
  return out;
}
