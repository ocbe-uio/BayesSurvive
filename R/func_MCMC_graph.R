#' @title Function to learn MRF graph
#'
#' @description
#' This an internal function for MCMC sampling
#'
#' @name func_MCMC_graph
#'
#' @import stats
#'
#' @param sobj a list containing observed data from \code{n} subjects;
#' \code{t}, \code{di}, \code{X}. See details for more information
#' @param hyperpar a list containing prior parameter values
#' @param ini a list containing prior parameters' ini values
#' @param S the number of subgroups
#' @param method a method option from
#' \code{c("Pooled", "CoxBVSSL", "Sub-struct")}
#' @param MRF_2b two different b in MRF prior for subgraphs G_ss and G_rs
#'
#' @return A list object with components "Sig" the updated covariance matrices,
#' "G.ini" the updated graph, "V.ini" the updated variances for precision
#' matrices in all subgroups, "C.ini" the updated precision matrices omega for
#' each subgroup
#'
#'
#' @export
func_MCMC_graph <- function(sobj, hyperpar, ini, S, method, MRF_2b) {
  n <- sobj$n
  p <- sobj$p
  SSig <- sobj$SSig

  pii <- hyperpar$pi.G # prior probability of edge inclusion
  a <- hyperpar$a # hyperparameter in MRF prior
  b <- hyperpar$b # hyperparameter in MRF prior
  lambda <- hyperpar$lambda # hyperparameter in precision matrix prior
  V0 <- hyperpar$V0 # (p x p) matrix of small variances in prior of precision matrix for each subgroup
  V1 <- hyperpar$V1 # (p x p) matrix of large variances in prior of precision matrix for each subgroup

  # Parameters to update:
  G <- G.MRF <- ini$G.ini # graph, large (pS x pS) matrix
  V <- ini$V.ini # list of length S with (p x p) matrices of updated variances for precision matrices in all subgroups
  C <- ini$C.ini # list of length S with (p x p) precision matrices omega for each subgroup
  Sig <- ini$Sig.ini # list of length S with (p x p) covariance (=inverse precision) matrices sigma

  gamma.ini <- ini$gamma.ini # list of length S with (p x 1) vectors of variable inclusion indicators in each subgroup
  gamma.ini <- unlist(gamma.ini) # write all gamma_s,i in one long vector

  # two different values for b in MRF prior for subgraphs G_ss and G_rs
  if (MRF_2b) {
    for (g in 1:S) { # b1 * G_ss
      G.MRF[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <-
        b[1] * G.MRF[(g - 1) * p + (1:p), (g - 1) * p + (1:p)]
    }
    for (g in 1:(S - 1)) { # b2 * G_rs
      for (r in g:(S - 1)) {
        G.MRF[(g - 1) * p + (1:p), r * p + (1:p)] <-
          G.MRF[r * p + (1:p), (g - 1) * p + (1:p)] <-
          b[2] * G.MRF[r * p + (1:p), (g - 1) * p + (1:p)]
      }
    }
  } else { # one value for b in MRF prior for all subgraphs
    G.MRF <- G.MRF * b
    b <- rep(b, 2)
  }

  # Update of precision matrix and graph within each subgroup
  # (analogous to SSSL algorithm for concentration (precision) graph models in function 'BayesGGM_SSVS_FixedV0V1' (Wang, 2015))

  for (g in 1:S) { # loop through subgroups

    V_g <- V[[g]] # variances for elements in precision matrix Omega_gg
    C_g <- C[[g]] # precision matrix Omega_gg in subgroup g
    G_g <- G[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] # subgraph G_gg in subgroup g
    n_g <- n[[g]] # sample size in subgroup g
    S_g <- SSig[[g]] # sample covariance matrix S=X'X in subgroup g
    Sig_g <- Sig[[g]] # covariance matrix Sigma_gg in subgroup g

    for (i in 1:p) { # loop through genes

      ind_noi <- setdiff(1:p, i)
      v_temp <- V_g[ind_noi, i]

      Sig11 <- Sig_g[ind_noi, ind_noi]
      Sig12 <- Sig_g[ind_noi, i]

      invC11 <- Sig11 - Sig12 %*% t(Sig12) / Sig_g[i, i] # Omega_11^(-1)

      Ci <- (S_g[i, i] + lambda) * invC11 + diag(1 / v_temp) # C^(-1)

      Ci <- (Ci + t(Ci)) / 2
      Ci_chol <- chol(Ci)

      mu_i <- -solve(Ci_chol, solve(t(Ci_chol), S_g[ind_noi, i]))
      beta <- mu_i + solve(Ci_chol, rnorm(p - 1))

      # Update of last column in Omega_gg
      C_g[ind_noi, i] <- C_g[i, ind_noi] <- beta # omega_12

      a_gam <- 0.5 * n_g + 1
      b_gam <- (S_g[i, i] + lambda) * 0.5
      gam <- rgamma(1, shape = a_gam, scale = 1 / b_gam)

      c <- t(beta) %*% invC11 %*% beta
      C_g[i, i] <- gam + c # omega_22

      # Below updating covariance matrix according to one-column change of precision matrix
      invC11beta <- invC11 %*% beta

      Sig_g[ind_noi, ind_noi] <- invC11 + invC11beta %*% t(invC11beta) / gam
      Sig_g[ind_noi, i] <- Sig_g[i, ind_noi] <- -invC11beta / gam
      Sig_g[i, i] <- 1 / gam

      # Update of variance matrix V_g and subgraph G_g

      if (i < p) {
        beta2 <- numeric(p)
        beta2[ind_noi] <- beta

        for (j in (i + 1):p) {
          # G where g_ss,ij=1 or 0:
          G.prop1 <- G.prop0 <- G.MRF
          G.prop1[(g - 1) * p + j, (g - 1) * p + i] <-
            G.prop1[(g - 1) * p + i, (g - 1) * p + j] <- b[1] # 1
          G.prop0[(g - 1) * p + j, (g - 1) * p + i] <-
            G.prop0[(g - 1) * p + i, (g - 1) * p + j] <- 0

          v0 <- V0[j, i]
          v1 <- V1[j, i]

          w1 <- -0.5 * log(v1) - 0.5 * beta2[j]^2 / v1 + log(pii)
          w2 <- -0.5 * log(v0) - 0.5 * beta2[j]^2 / v0 + log(1 - pii)

          wa <- w1 + (a * sum(gamma.ini) + t(gamma.ini) %*% G.prop1 %*% gamma.ini)
          wb <- w2 + (a * sum(gamma.ini) + t(gamma.ini) %*% G.prop0 %*% gamma.ini)

          w_max <- max(wa, wb)

          w <- exp(wa - w_max) / (exp(wa - w_max) + exp(wb - w_max))

          z <- (runif(1) < as.numeric(w)) # acceptance/ rejection of proposal
          v <- ifelse(z, v1, v0)
          V_g[j, i] <- V_g[i, j] <- v

          G_g[j, i] <- G_g[i, j] <- as.numeric(z)
        }
      }
    }

    V[[g]] <- V_g
    C[[g]] <- C_g
    G[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <- G_g
    Sig[[g]] <- Sig_g
  }

  ####

  # Update of the subgraph between subgroups (only for CoxBVS-SL model; for Sub-struct this subgraph
  # remains fixed to starting values)

  if (method == "CoxBVSSL") {
    for (g in 1:(S - 1)) { # loop through subgroups
      for (r in g:(S - 1)) { # loop through subgroups

        G_rs <- G[(g - 1) * p + (1:p), r * p + (1:p)] # subgraph between subgroups g and r>g

        for (i in 1:p) { # loop through genes
          G.prop1 <- G.prop0 <- G.MRF # G where g_rs,ii=1 or 0
          G.prop1[(g - 1) * p + i, r * p + i] <- G.prop1[r * p + i, (g - 1) * p + i] <- b[2]
          G.prop0[(g - 1) * p + i, r * p + i] <- G.prop0[r * p + i, (g - 1) * p + i] <- 0


          # compute conditional distribution for proposal g_rs,ii=1
          wa <- log(pii) + (a * sum(gamma.ini) + t(gamma.ini) %*% G.prop1 %*% gamma.ini)
          wb <- log(1 - pii) + (a * sum(gamma.ini) + t(gamma.ini) %*% G.prop0 %*% gamma.ini)

          w_max <- max(wa, wb)
          pg <- exp(wa - w_max) / (exp(wa - w_max) + exp(wb - w_max))

          G_rs[i, i] <- as.numeric(runif(1) < pg) # acceptance/ rejection of proposal
        }
        G[(g - 1) * p + (1:p), r * p + (1:p)] <-
          G[r * p + (1:p), (g - 1) * p + (1:p)] <- G_rs # update of subgraph
      }
    }
  }

  return(list(Sig.ini = Sig, G.ini = G, V.ini = V, C.ini = C))
}
