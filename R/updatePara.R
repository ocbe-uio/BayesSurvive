#' @title Subfunctions to update parameters
#'
#' @description
#' This contains subfunctions to update parameters gammas, betas, baseline
#' hazard and graph learning parameters
#'
#' @name UpdateGamma
#'
#' @importFrom mvtnorm dmvnorm
#'
#' @param sobj a list containing observed data
#' @param hyperpar a list containing prior parameter values
#' @param ini a list containing prior parameters' initial values
#' @param S the number of subgroups
#' @param method a method option from
#' \code{c("Pooled", "CoxBVSSL", "Sub-struct", "Subgroup")}
#' @param MRF_G logical value. \code{MRF_G = TRUE} is to fix the MRF graph which
#' is provided in the argument \code{hyperpar}, and \code{MRF_G = FALSE} is to
#' use graphical model for learning the MRF graph
#' @param MRF_2b two different b in MRF prior for subgraphs G_ss and G_rs
#' @inheritParams func_MCMC
#'
#' @return A list object with two components for the latent variable selection
#' indicators gamma with either independent Bernoulli prior
# (standard approaches) or with MRF prior
#'
#' @export
UpdateGamma <- function(sobj, hyperpar, ini, S, method, MRF_G, MRF_2b, cpp = FALSE) {
  # Update latent variable selection indicators gamma with either independent Bernoulli prior
  # (standard approaches) or with MRF prior.
  # browser()
  if (cpp) {
    return(UpdateGamma_cpp(sobj, hyperpar, ini, S, method, MRF_G, MRF_2b))
  }
  p <- sobj$p
  tau <- hyperpar$tau
  cb <- hyperpar$cb
  pi <- hyperpar$pi.ga
  a <- hyperpar$a
  b <- hyperpar$b

  beta.ini <- ini$beta.ini
  gamma.ini <- ini$gamma.ini

  if (method %in% c("Pooled") && MRF_G) {
    G.ini <- hyperpar$G
  }

  # if (method %in% c("CoxBVSSL", "Sub-struct") ||
  #     (method == "Pooled" && !MRF_G)) {
  if (!MRF_G) {
    G.ini <- ini$G.ini
  }

  # two different b in MRF prior for subgraphs G_ss and G_rs
  if (MRF_2b && !MRF_G) {
    for (g in 1:S) { # b1*G_ss
      G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <-
        b[1] * G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)]
    }
    for (g in 1:(S - 1)) { # b2*G_rs
      for (r in g:(S - 1)) {
        G.ini[(g - 1) * p + (1:p), r * p + (1:p)] <-
          G.ini[r * p + (1:p), (g - 1) * p + (1:p)] <-
          b[2] * G.ini[r * p + (1:p), (g - 1) * p + (1:p)]
      }
    }
  } else {
    # if (method %in% c("CoxBVSSL", "Sub-struct") ||
    #     (method == "Pooled" && !MRF_G)) {
    if (!MRF_G) {
      G.ini <- G.ini * b
    }
  }

  if (method == "Pooled" && MRF_G) {
    post.gamma <- rep(0, p)

    for (j in 1:p) {
      # wa   = dnorm(beta.ini[j], mean = 0, sd = cb*tau) * pi
      # wb   = dnorm(beta.ini[j], mean = 0, sd = tau) * (1 - pi)
      # pgam = wa/(wa + wb)
      # u = runif(1)
      # gamma.ini[j]  = ifelse(u < pgam, 1, 0)
      # post.gamma[j] = pgam



      beta <- beta.ini[j]

      ga.prop1 <- ga.prop0 <- gamma.ini # gamma with gamma_g,j=1 or 0
      ga.prop1[j] <- 1
      ga.prop0[j] <- 0
      ga.prop1 <- unlist(ga.prop1)
      ga.prop0 <- unlist(ga.prop0)

      wa <- (a * sum(ga.prop1) + t(ga.prop1) %*% G.ini %*% ga.prop1) +
        dnorm(beta, mean = 0, sd = tau * cb, log = TRUE)
      wb <- (a * sum(ga.prop0) + t(ga.prop0) %*% G.ini %*% ga.prop0) +
        dnorm(beta, mean = 0, sd = tau, log = TRUE)

      w_max <- max(wa, wb)
      pg <- exp(wa - w_max) / (exp(wa - w_max) + exp(wb - w_max))

      gamma.ini[j] <- as.numeric(runif(1) < pg)
      post.gamma[j] <- pg
    }
  } else {
    post.gamma <- rep(list(rep(0, p)), S)

    ## ? "Subgroup" might be not needed
    # if (method == "Subgroup") {
    if (MRF_G) {
      for (g in 1:S) { # loop through subgroups
        for (j in 1:p) {
          wa <- dnorm((beta.ini[[g]])[j], mean = 0, sd = cb * tau) * pi
          wb <- dnorm((beta.ini[[g]])[j], mean = 0, sd = tau) * (1 - pi)
          pgam <- wa / (wa + wb)
          u <- runif(1)
          gamma.ini[[g]][j] <- ifelse(u < pgam, 1, 0)
          post.gamma[[g]][j] <- pgam
        }
      }
    } else { # CoxBVS-SL or Sub-struct model

      for (g in 1:S) { # loop through subgroups
        for (j in 1:p) {
          beta <- (beta.ini[[g]])[j]

          ga.prop1 <- ga.prop0 <- gamma.ini # gamma with gamma_g,j=1 or 0
          ga.prop1[[g]][j] <- 1
          ga.prop0[[g]][j] <- 0
          ga.prop1 <- unlist(ga.prop1)
          ga.prop0 <- unlist(ga.prop0)

          wa <- (a * sum(ga.prop1) + t(ga.prop1) %*% G.ini %*% ga.prop1) +
            dnorm(beta, mean = 0, sd = tau * cb, log = TRUE)
          wb <- (a * sum(ga.prop0) + t(ga.prop0) %*% G.ini %*% ga.prop0) +
            dnorm(beta, mean = 0, sd = tau, log = TRUE)

          w_max <- max(wa, wb)
          pg <- exp(wa - w_max) / (exp(wa - w_max) + exp(wb - w_max))

          gamma.ini[[g]][j] <- as.numeric(runif(1) < pg)
          post.gamma[[g]][j] <- pg
        }
      }
    }
  }

  return(list(gamma.ini = gamma.ini, post.gamma = post.gamma))
}
# the end of "UpdateGamma" function

####

# Update joint posterior distribution

calJpost <- function(sobj, hyperpar, ini, S, method, MRF_G, MRF_2b, cpp) {
  if (cpp && method == "Pooled" && MRF_G) { # TODO: keep only cpp once translation is complete
    return(calJpost_cpp(sobj, hyperpar, ini, S, method, MRF_G, MRF_2b))
  }
  # hyperparameters
  p <- sobj$p
  c0 <- hyperpar$c0
  pi.ga <- hyperpar$pi.ga
  tau <- hyperpar$tau
  cb <- hyperpar$cb

  if (method %in% c("CoxBVSSL", "Sub-struct") ||
    (method == "Pooled" && !MRF_G)) {
    lambda <- hyperpar$lambda
    a <- hyperpar$a
    b <- hyperpar$b
    pi.G <- hyperpar$pi.G
    G.ini <- ini$G.ini
  }

  if (method == "Pooled" && MRF_G) {
    n <- sobj$n
    X <- sobj$X
    J <- hyperpar$J
    ind.r_d <- hyperpar$ind.r_d
    ind.d <- hyperpar$ind.d
    hPriorSh <- hyperpar$hPriorSh

    gamma.ini <- ini$gamma.ini
    beta.ini <- ini$beta.ini
    h <- ini$h
    cbtau <- tau * ifelse(gamma.ini == 1, cb, 1)

    # erg <- calJpost.helper(cbtau, X, beta.ini, h, hPriorSh, c0, n, p, J, ind.r_d, ind.d)
    erg <- calJpost_helper_cpp(cbtau, X, beta.ini, h, hPriorSh, c0, ind.r_d, ind.d)
    loglike <- erg$loglike1
    logpriorBeta <- erg$logpriorBeta1
    logpriorH <- erg$logpriorH1

    logpriorGamma <- sum(gamma.ini * log(pi.ga)) + sum((1 - gamma.ini) * log(1 - pi.ga))
    logjpost <- loglike + logpriorGamma + logpriorBeta + logpriorH
  } else {
    loglike <- logpriorBeta <- logpriorH <- logpriorGamma <-
      logjpost <- logpriorOmega <- logpriorX <- numeric()

    for (g in 1:S) {
      n <- sobj$n[[g]]
      X <- sobj$X[[g]]
      J <- hyperpar$J[[g]]
      ind.r_d <- hyperpar$ind.r_d[[g]]
      ind.d <- hyperpar$ind.d[[g]]
      hPriorSh <- hyperpar$hPriorSh[[g]]

      gamma.ini <- ini$gamma.ini[[g]]
      beta.ini <- ini$beta.ini[[g]]
      h <- ini$h[[g]]
      cbtau <- tau * ifelse(gamma.ini == 1, cb, 1)

      # erg <- calJpost.helper(cbtau, X, beta.ini, h, hPriorSh, c0, n, p, J, ind.r_d, ind.d)
      erg <- calJpost_helper_cpp(cbtau, X, beta.ini, h, hPriorSh, c0, ind.r_d, ind.d)
      loglike[g] <- erg$loglike1
      logpriorBeta[g] <- erg$logpriorBeta1
      logpriorH[g] <- erg$logpriorH1

      # if (method == "Subgroup") {
      if (MRF_G) {
        logpriorGamma[g] <- sum(gamma.ini * log(pi.ga)) +
          sum((1 - gamma.ini) * log(1 - pi.ga))
        logjpost[g] <- loglike[g] + logpriorGamma[g] + logpriorBeta[g] + logpriorH[g]
      } else { # CoxBVSSL/ Sub-struct model

        C.ini <- ini$C.ini[[g]]
        V.ini <- ini$V.ini[[g]]
        Sig.ini <- ini$Sig.ini[[g]]

        omega.mat <- matrix(0, p, p)
        for (i in 1:p) {
          omega.mat[i, ] <- dnorm(C.ini[i, ],
            mean = rep(0, p),
            sd = sqrt(V.ini[i, ]), log = TRUE
          )
        }
        diag(omega.mat) <- dexp(diag(C.ini), rate = lambda / 2, log = TRUE)

        logpriorOmega[g] <- sum(omega.mat[upper.tri(omega.mat, diag = TRUE)])
        logpriorX[g] <- sum(dmvnorm(X,
          mean = rep(0, p),
          sigma = Sig.ini, log = TRUE
        ))
      }
    }
  }
  # if (method %in% c("CoxBVSSL", "Sub-struct") ||
  #     (method == "Pooled" && !MRF_G)) {
  if (!MRF_G) {
    pii.mat <- matrix(0, p * S, p * S)

    # id.mat is "full" graph (with all possible edges set to 1)
    if (method == "CoxBVSSL") {
      id.mat <- do.call("cbind", rep(list(
        do.call("rbind", rep(list(diag(1, p, p)), S))
      ), S))
    } else {
      id.mat <- matrix(0, p * S, p * S)
    }
    for (g in 1:S) {
      id.mat[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <- matrix(1, p, p)
    }

    pii.mat[id.mat == 1 & G.ini == 0] <- log(1 - pi.G)
    pii.mat[id.mat == 1 & G.ini == 1] <- log(pi.G)

    logpriorGraph <- sum(pii.mat[upper.tri(pii.mat, diag = FALSE)])

    gamma.vec <- unlist(ini$gamma.ini)

    if (MRF_2b) {
      for (g in 1:S) { # b1 * G_ss
        G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <-
          b[1] * G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)]
      }
      for (g in 1:(S - 1)) { # b2 * G_rs
        for (r in g:(S - 1)) {
          G.ini[(g - 1) * p + (1:p), r * p + (1:p)] <-
            G.ini[r * p + (1:p), (g - 1) * p + (1:p)] <-
            b[2] * G.ini[r * p + (1:p), (g - 1) * p + (1:p)]
        }
      }
      b <- 1
    }
    logpriorGamma <- (a * sum(gamma.vec) + b * t(gamma.vec) %*% G.ini %*% gamma.vec)

    logjpost <- sum(loglike) + sum(logpriorBeta) + sum(logpriorH) +
      logpriorGamma + sum(logpriorOmega) + sum(logpriorX) + logpriorGraph
  }
  return(list(loglike = loglike, logjpost = logjpost))
}
