#' BayesSurvive
#' @title Subfunctions to update parameters
#'
#' @description
#' This contains subfunctions to update parameters gammas, betas, baseline
#' hazard and graph learning parameters
#'
#' @name UpdateGamma
#'
#' @param sobj a list containing observed data
#' @param hyperpar a list containing prior parameter values
#' @param ini a list containing prior parameters' initial values
#' @param S the number of subgroups
#' @param method a method option from
#' \code{c("Pooled", "CoxBVSSL", "Sub-struct", "Subgroup")}
#' @param MRF_G logical value. \code{MRF_G = TRUE} is to fix the MRF graph which
#' is provided in the argument \code{hyperpar}, and \code{MRF_G = FALSE} is to
#' use graphical model for leanring the MRF graph
#' @param MRF_2b two different b in MRF prior for subgraphs G_ss and G_rs
#'
#' @return An object of ...
#'
#' @examples
#'
#' # Load the example dataset
#'
#' @export
UpdateGamma <- function(sobj, hyperpar, ini, S, method, MRF_G, MRF_2b) {
  # Update latent variable selection indicators gamma with either independent Bernoulli prior
  # (standard approaches) or with MRF prior.

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
      G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <- b[1] * G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)]
    }
    for (g in 1:(S - 1)) { # b2*G_rs
      for (r in g:(S - 1)) {
        G.ini[(g - 1) * p + (1:p), r * p + (1:p)] <- G.ini[r * p + (1:p), (g - 1) * p + (1:p)] <- b[2] * G.ini[r * p + (1:p), (g - 1) * p + (1:p)]
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

      wa <- (a * sum(ga.prop1) + t(ga.prop1) %*% G.ini %*% ga.prop1) + dnorm(beta, mean = 0, sd = tau * cb, log = TRUE)
      wb <- (a * sum(ga.prop0) + t(ga.prop0) %*% G.ini %*% ga.prop0) + dnorm(beta, mean = 0, sd = tau, log = TRUE)

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

          wa <- (a * sum(ga.prop1) + t(ga.prop1) %*% G.ini %*% ga.prop1) + dnorm(beta, mean = 0, sd = tau * cb, log = TRUE)
          wb <- (a * sum(ga.prop0) + t(ga.prop0) %*% G.ini %*% ga.prop0) + dnorm(beta, mean = 0, sd = tau, log = TRUE)

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

calJpost <- function(sobj, hyperpar, ini, S, method, MRF_G, MRF_2b) {
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
    erg <- calJpost_helper_cpp(cbtau, X, beta.ini, h, hPriorSh, c0, J, ind.r_d, ind.d)
    loglike <- erg$loglike1
    logpriorBeta <- erg$logpriorBeta1
    logpriorH <- erg$logpriorH1

    logpriorGamma <- sum(gamma.ini * log(pi.ga)) + sum((1 - gamma.ini) * log(1 - pi.ga))
    logjpost <- loglike + logpriorGamma + logpriorBeta + logpriorH
  } else {
    loglike <- logpriorBeta <- logpriorH <- logpriorGamma <- logjpost <- logpriorOmega <- logpriorX <- numeric()

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
      erg <- calJpost_helper_cpp(cbtau, X, beta.ini, h, hPriorSh, c0, J, ind.r_d, ind.d)
      loglike[g] <- erg$loglike1
      logpriorBeta[g] <- erg$logpriorBeta1
      logpriorH[g] <- erg$logpriorH1

      # if (method == "Subgroup") {
      if (MRF_G) {
        logpriorGamma[g] <- sum(gamma.ini * log(pi.ga)) + sum((1 - gamma.ini) * log(1 - pi.ga))
        logjpost[g] <- loglike[g] + logpriorGamma[g] + logpriorBeta[g] + logpriorH[g]
      } else { # CoxBVSSL/ Sub-struct model

        C.ini <- ini$C.ini[[g]]
        V.ini <- ini$V.ini[[g]]
        Sig.ini <- ini$Sig.ini[[g]]

        omega.mat <- matrix(0, p, p)
        for (i in 1:p) {
          omega.mat[i, ] <- dnorm(C.ini[i, ], mean = rep(0, p), sd = sqrt(V.ini[i, ]), log = TRUE)
        }
        diag(omega.mat) <- dexp(diag(C.ini), rate = lambda / 2, log = TRUE)

        logpriorOmega[g] <- sum(omega.mat[upper.tri(omega.mat, diag = TRUE)])
        logpriorX[g] <- sum(dmvnorm(X, mean = rep(0, p), sigma = Sig.ini, log = TRUE))
      }
    }
  }
  # if (method %in% c("CoxBVSSL", "Sub-struct") ||
  #     (method == "Pooled" && !MRF_G)) {
  if (!MRF_G) {
    pii.mat <- matrix(0, p * S, p * S)

    # id.mat is "full" graph (with all possible edges set to 1)
    if (method == "CoxBVSSL") {
      id.mat <- do.call("cbind", rep(list(do.call("rbind", rep(list(diag(1, p, p)), S))), S))
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
        G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <- b[1] * G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)]
      }
      for (g in 1:(S - 1)) { # b2 * G_rs
        for (r in g:(S - 1)) {
          G.ini[(g - 1) * p + (1:p), r * p + (1:p)] <- G.ini[r * p + (1:p), (g - 1) * p + (1:p)] <- b[2] * G.ini[r * p + (1:p), (g - 1) * p + (1:p)]
        }
      }
      b <- 1
    }
    logpriorGamma <- (a * sum(gamma.vec) + b * t(gamma.vec) %*% G.ini %*% gamma.vec)

    logjpost <- sum(loglike) + sum(logpriorBeta) + sum(logpriorH) + logpriorGamma + sum(logpriorOmega) + sum(logpriorX) + logpriorGraph
  }
  return(list(loglike = loglike, logjpost = logjpost))
}


####

# # Helper function for "UpdateRP.lee11" function.
#
# UpdateRP.lee11.helper <- function(n, p, x, J, ind.r, ind.d, ind.r_d, be.ini,
#                                   ga.ini, h, tau, cb) {
#   acceptl <- rep(0, p)
#
#   for (j in 1:p) {
#     sd.be <- tau * ifelse(ga.ini[j] == 1, cb, 1)
#
#     xbeta <- apply(x * matrix(rep(be.ini, n), n, p, byrow = T), 1, sum)
#     xbeta[xbeta > 700] <- 700
#     exp.xbeta <- exp(xbeta)
#     x.exp.xbeta <- x[, j] * exp.xbeta
#     exp.xbeta.mat <- matrix(rep(exp.xbeta, J), n, J)
#     x.exp.xbeta.mat <- matrix(rep(x.exp.xbeta, J), n, J)
#     D1.1st <- -h * colSums(x.exp.xbeta.mat * ind.r_d)
#
#     h.mat <- matrix(rep(h, n), n, J, byrow = T)
#     h.exp.xbeta.mat <- -h.mat * exp.xbeta.mat
#     h.exp.xbeta.mat[h.exp.xbeta.mat > -10^(-7)] <- -10^(-7)
#     exp.h.exp.xbeta.mat <- exp(h.exp.xbeta.mat)
#
#     D1.2nd.den <- 1 - exp.h.exp.xbeta.mat
#     D1.2nd.num <- exp.h.exp.xbeta.mat * x.exp.xbeta.mat
#     D1.2nd <- h * colSums(D1.2nd.num / D1.2nd.den * ind.d)
#     D1 <- sum(D1.1st + D1.2nd) - 1 / sd.be^2 * be.ini[j]
#
#     x.sq.exp.xbeta <- x[, j]^2 * exp.xbeta
#     x.sq.exp.xbeta.mat <- matrix(rep(x.sq.exp.xbeta, J), n, J)
#
#     D2.1st <- -h * colSums(x.sq.exp.xbeta.mat * ind.r_d)
#     D2.2nd.den <- D1.2nd.den^2
#     D2.2nd.num <- exp.h.exp.xbeta.mat * x.sq.exp.xbeta.mat * (1 - exp.h.exp.xbeta.mat + h.exp.xbeta.mat)
#     D2.2nd <- h * colSums(D2.2nd.num / D2.2nd.den * ind.d)
#     D2 <- sum(D2.1st + D2.2nd) - 1 / sd.be^2
#
#     be.prop.me <- be.ini[j] - D1 / D2
#     be.prop.var <- -2.4^2 / D2
#     be.prop <- be.ini
#     be.prop[j] <- rnorm(1, mean = be.prop.me, sd = sqrt(be.prop.var))
#
#     xbeta.prop <- xbeta - x[, j] * be.ini[j] + x[, j] * be.prop[j]
#     xbeta.prop[xbeta.prop > 700] <- 700
#     exp.xbeta.prop <- exp(xbeta.prop)
#     x.exp.xbeta.prop <- x[, j] * exp.xbeta.prop
#     exp.xbeta.prop.mat <- matrix(rep(exp.xbeta.prop, J), n, J)
#     x.exp.xbeta.prop.mat <- matrix(rep(x.exp.xbeta.prop, J), n, J)
#     D1.1st.prop <- -h * colSums(x.exp.xbeta.prop.mat * ind.r_d)
#
#     h.exp.xbeta.prop.mat <- -h.mat * exp.xbeta.prop.mat
#     h.exp.xbeta.prop.mat[h.exp.xbeta.prop.mat > -10^(-7)] <- -10^(-7)
#     exp.h.exp.xbeta.prop.mat <- exp(h.exp.xbeta.prop.mat)
#
#     D1.2nd.den.prop <- 1 - exp.h.exp.xbeta.prop.mat
#     D1.2nd.num.prop <- exp.h.exp.xbeta.prop.mat * x.exp.xbeta.prop.mat
#     D1.2nd.prop <- h * colSums(D1.2nd.num.prop / D1.2nd.den.prop * ind.d)
#     D1.prop <- sum(D1.1st.prop + D1.2nd.prop) - 1 / sd.be^2 * be.prop[j]
#
#     x.sq.exp.xbeta.prop <- x[, j]^2 * exp.xbeta.prop
#     x.sq.exp.xbeta.prop.mat <- matrix(rep(x.sq.exp.xbeta.prop, J), n, J)
#
#     D2.1st.prop <- -h * colSums(x.sq.exp.xbeta.prop.mat * ind.r_d)
#     D2.2nd.den.prop <- D1.2nd.den.prop^2
#     D2.2nd.num.prop <- exp.h.exp.xbeta.prop.mat * x.sq.exp.xbeta.prop.mat * (1 - exp.h.exp.xbeta.prop.mat + h.exp.xbeta.prop.mat)
#     D2.2nd.prop <- h * colSums(D2.2nd.num.prop / D2.2nd.den.prop * ind.d)
#     D2.prop <- sum(D2.1st.prop + D2.2nd.prop) - 1 / sd.be^2
#
#     be.prop.me.ini <- be.prop[j] - D1.prop / D2.prop
#     be.prop.var.ini <- -2.4^2 / D2.prop
#
#     first.sum <- colSums(exp.xbeta.mat * ind.r_d)
#     second.sum <- colSums(log(D1.2nd.den) * ind.d)
#     loglh.ini <- sum(-h * first.sum + second.sum)
#
#     first.sum.prop <- colSums(exp.xbeta.prop.mat * ind.r_d)
#     second.sum.prop <- colSums(log(D1.2nd.den.prop) * ind.d)
#     loglh.prop <- sum(-h * first.sum.prop + second.sum.prop)
#
#     logprior.prop <- dnorm(be.prop[j], mean = 0, sd = sd.be, log = TRUE)
#     logprior.ini <- dnorm(be.ini[j], mean = 0, sd = sd.be, log = TRUE)
#     logprop.prop <- dnorm(be.prop[j], mean = be.prop.me.ini, sd = sqrt(be.prop.var.ini), log = TRUE)
#     logprop.ini <- dnorm(be.ini[j], mean = be.prop.me, sd = sqrt(be.prop.var), log = TRUE)
#
#     logR <- loglh.prop - loglh.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop
#     u <- log(runif(1)) < logR
#     if (u == 1) {
#       be.ini[j] <- be.prop[j]
#       acceptl[j] <- acceptl[j] + 1
#     }
#   }
#   return(list("be.ini" = be.ini, "acceptl" = acceptl))
# }
# # the end of "UpdateRP.lee11.helper" function

# # Update regression coefficients according to Lee et al. (2011)
#
# UpdateRP.lee11 <- function(sobj, hyperpar, ini, S, method)
# {
#   p   = sobj$p
#   tau = hyperpar$tau
#   cb  = hyperpar$cb
#
#   if(method == "Pooled"){
#     n       = sobj$n
#     x       = sobj$X
#     J       = hyperpar$J
#     ind.r   = hyperpar$ind.r
#     ind.d   = hyperpar$ind.d
#     ind.r_d = hyperpar$ind.r_d
#     be.ini  = ini$beta.ini
#     ga.ini  = ini$gamma.ini
#     h       = ini$h
#
#     #erg = UpdateRP.lee11.helper(n, p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
#     erg = updateRP_genomic_cpp(p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
#
#     beta.ini  = erg$be.ini
#     acceptlee = erg$acceptl
#
#   }else{
#     beta.ini = acceptlee = vector("list", S )
#     for(g in 1:S){  # loop through subgroups
#
#       n       = sobj$n[[g]]
#       x       = sobj$X[[g]]
#       J       = hyperpar$J[[g]]
#       ind.r   = hyperpar$ind.r[[g]]
#       ind.d   = hyperpar$ind.d[[g]]
#       ind.r_d = hyperpar$ind.r_d[[g]]
#       be.ini  = ini$beta.ini[[g]]
#       ga.ini  = ini$gamma.ini[[g]]
#       h       = ini$h[[g]]
#
#       #erg = UpdateRP.lee11.helper(n, p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
#       erg = updateRP_genomic_cpp(p, x, J, ind.r, ind.d, ind.r_d, be.ini, ga.ini, h, tau, cb)
#
#       beta.ini[[g]] = erg$be.ini
#       acceptlee[[g]] = erg$acceptl
#     }
#   }
#   return(list(beta.ini = beta.ini, acceptlee = acceptlee))
# }
# # the end of "UpdateRP.lee11" function

####

# # Helper function for "calJpost" function.
#
# calJpost.helper <- function(cbtau, X, beta.ini, h, hPriorSh, c0, n, p, J, ind.r_d, ind.d) {
#   xbeta <- apply(X * matrix(rep(beta.ini, n), n, p, byrow = T), 1, sum)
#   xbeta.mat <- matrix(rep(xbeta, J), n, J)
#   xbeta.dum <- xbeta
#   xbeta.dum[xbeta.dum > 700] <- 700
#   exp.xbeta <- exp(xbeta.dum)
#   exp.xbeta.mat <- matrix(rep(exp.xbeta, J), n, J)
#   h.mat <- matrix(rep(h, n), n, J, byrow = T)
#   h.exp.xbeta.mat <- -h.mat * exp.xbeta.mat
#   h.exp.xbeta.mat[h.exp.xbeta.mat > -10^(-7)] <- -10^(-7)
#
#   first.sum.ini <- sum(-h * apply(exp.xbeta.mat * ind.r_d, 2, sum))
#   second.sum.ini <- sum(apply(log(1 - exp(h.exp.xbeta.mat)) * ind.d, 2, sum))
#   loglike1 <- first.sum.ini + second.sum.ini
#   logpriorBeta1 <- sum(dnorm(beta.ini, mean = rep(0, p), sd = cbtau, log = TRUE))
#   logpriorH1 <- sum(dgamma(h, shape = hPriorSh, rate = c0, log = TRUE))
#
#   return(list("loglike1" = loglike1, "logpriorBeta1" = logpriorBeta1, "logpriorH1" = logpriorH1))
# }


# setting.interval = function(y, delta, s, J)
# {
#   ind.d = ind.r = matrix(0, length(y), J)
#
#   for(i in which(delta == 1)){
#     d.mat.ind = min(which(s - y[i] >=0))
#     ind.d[i, d.mat.ind] = 1
#     ind.r[i, 1:d.mat.ind] = 1
#   }
#
#   for(i in which(delta == 0)){
#     r.mat.ind = min(which(s - y[i] >=0))
#     ind.r[i, 1:r.mat.ind] = 1
#   }
#
#   ind.r[,1] = 1
#   ind.r_d	= ind.r - ind.d;
#
#   d = apply(ind.d, 2, sum)
#
#   return(list(s = s, J = J, ind.r = ind.r, ind.d = ind.d, d = d, ind.r_d = ind.r_d))
# }
