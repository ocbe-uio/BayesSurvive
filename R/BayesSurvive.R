#' @title Fit Bayesian Cox Models
#'
#' @description
#' This is the main function to fit a Bayesian Cox model with graph-structured
#' selection priors for sparse identification of high-dimensional covariates.
#'
#' @name BayesSurvive
#'
#' @importFrom Rcpp evalCpp
#' @importFrom survival survreg
#' @importFrom stats runif
#' @importFrom methods is
#'
#' @param survObj a list containing observed data from \code{n} subjects with
#' components \code{t}, \code{di}, \code{X}. For graphical learning of the
#' Markov random field prior, \code{survObj} should be a list of the list with
#' survival and covariates data. For subgroup models with or without graphical
#' learning, \code{survObj} should be a list of multiple lists with each
#' component list representing each subgroup's survival and covariates data
#' @param model.type a method option from
#' \code{c("Pooled", "CoxBVSSL", "Sub-struct")}. To enable graphical learning
#' for "Pooled" model, please specify \code{list(survObj)} where \code{survObj}
#' is the list of \code{t}, \code{di} and \code{X}
#' @param MRF2b logical value. \code{MRF2b = TRUE} means two different
#' hyperparameters b in MRF prior (values b01 and b02) and \code{MRF2b = FALSE}
#' means one hyperparamter b in MRF prior
#' @param MRF.G logical value. \code{MRF.G = TRUE} is to fix the MRF graph which
#' is provided in the argument \code{hyperpar}, and \code{MRF.G = FALSE} is to
#' use graphical model for leanring the MRF graph
#' @param g.ini initial values for latent edge inclusion indicators in graph,
#' should be a value in [0,1]. 0 or 1: set all random edges to 0 or 1; value in
#' (0,1): rate of indicators randomly set to 1, the remaining indicators are 0
#' @param hyperpar a list containing prior parameter values
#' @param initial a list containing prior parameters' initial values
#' @param nIter the number of iterations of the chain
#' @param burnin number of iterations to discard at the start of the chain.
#' Default is 0
#' @param thin thinning MCMC intermediate results to be stored
#' @param output_graph_para allow (\code{TRUE}) or suppress (\code{FALSE}) the
#' output for parameters 'G', 'V', 'C' and 'Sig' in the graphical model
#' if \code{MRF.G = FALSE}
#' @param verbose logical value to display the progess of MCMC
#' @param cpp logical, whether to use C++ code for faster computation
#'
#'
#' @return An object of class \code{BayesSurvive} is saved as
#' \code{obj_BayesSurvive.rda} in the output file, including the following components:
#' \itemize{
#' \item input - a list of all input parameters by the user
#' \item output - a list of the all output estimates:
#' \itemize{
#' \item "\code{gamma.p}" - a matrix with MCMC intermediate estimates of the indicator variables of regression coefficients.
#' \item "\code{beta.p}" - a matrix with MCMC intermediate estimates of the regression coefficients.
#' \item "\code{h.p}" - a matrix with MCMC intermediate estimates of the increments in the cumulative baseline hazard in each interval.
#' }
#' \item call - the matched call.
#' }
#'
#'
#' @examples
#'
#' library("BayesSurvive")
#' set.seed(123)
#'
#' # Load the example dataset
#' data("simData", package = "BayesSurvive")
#'
#' dataset <- list(
#'   "X" = simData[[1]]$X,
#'   "t" = simData[[1]]$time,
#'   "di" = simData[[1]]$status
#' )
#'
#' # Initial value: null model without covariates
#' initial <- list("gamma.ini" = rep(0, ncol(dataset$X)))
#' # Hyperparameters
#' hyperparPooled <- list(
#'   "c0"     = 2, # prior of baseline hazard
#'   "tau"    = 0.0375, # sd (spike) for coefficient prior
#'   "cb"     = 20, # sd (spike) for coefficient prior
#'   "pi.ga"  = 0.02, # prior variable selection probability for standard Cox models
#'   "a"      = -4, # hyperparameter in MRF prior
#'   "b"      = 0.1, # hyperparameter in MRF prior
#'   "G"      = simData$G # hyperparameter in MRF prior
#' )
#'
#' \donttest{
#' # run Bayesian Cox with graph-structured priors
#' fit <- BayesSurvive(
#'   survObj = dataset, hyperpar = hyperparPooled,
#'   initial = initial, nIter = 100
#' )
#'
#' # show posterior mean of coefficients and 95% credible intervals
#' library("GGally")
#' plot(fit) +
#'   coord_flip() +
#'   theme(axis.text.x = element_text(angle = 90, size = 7))
#' }
#'
#' @export
BayesSurvive <- function(survObj,
                         model.type = "Pooled",
                         MRF2b = FALSE,
                         MRF.G = TRUE,
                         g.ini = 0,
                         hyperpar = NULL,
                         initial = NULL,
                         nIter = 1,
                         burnin = 0,
                         thin = 1,
                         output_graph_para = FALSE,
                         verbose = TRUE,
                         cpp = FALSE) {
  # Validation
  stopifnot(burnin < nIter)

  # same number of covariates p in all subgroups
  p <- ifelse(is.list(survObj[[1]]), NCOL(survObj[[1]]$X), NCOL(survObj$X))
  Beta.ini <- numeric(p)
  gamma.ini <- initial$gamma.ini
  if (sum(gamma.ini) != 0) {
    Beta.ini[which(gamma.ini == 0)] <- runif(sum(gamma.ini == 0), -0.02, 0.02)
    Beta.ini[which(gamma.ini == 1)] <- runif(sum(gamma.ini == 1), -2, 2)
  }
  beta.ini <- Beta.ini

  # for Pooled model combine/ pool data of all subgroups
  if (model.type == "Pooled" && MRF.G) {
    S <- 1 # number of subgroups
    survObj <- list(
      "t" = survObj$t, "di" = survObj$di, "X" = survObj$X,
      "SSig" = t(survObj$X) %*% survObj$X
    ) # sample covariance matrix
    survObj$n <- length(survObj$t)
    survObj$p <- p <- NCOL(survObj$X)

    fit <- survreg(Surv(survObj$t, survObj$di, type = c("right")) ~ 1,
      dist = "weibull", x = TRUE, y = TRUE
    )
    kappa0 <- 1 / exp(fit$icoef["Log(scale)"])
    eta0 <- exp(fit$coefficients)^(-kappa0)
    # Initial value: null model without covariates
    log.like <- coxph(Surv(survObj$t, survObj$di, type = c("right")) ~ 1)$loglik
  } else { # for Subgroup or CoxBVS-SL model keep the data of all subgroups separate

    S <- length(survObj) # number of subgroups
    for (g in 1:S) { # separate scaling of covariates
      survObj[[g]]$X <- scale(survObj[[g]]$X, scale = TRUE)
    }

    # estimate shape and scale parameter of Weibull distribution
    # (hyperparameters for prior of H* (mean function in Gamma process prior for baseline hazard))
    fit <- lapply(survObj, function(x) {
      survreg(Surv(x$t, x$di, type = c("right")) ~ 1, dist = "weibull", x = TRUE, y = TRUE)
    })
    kappa0 <- lapply(fit, function(x) {
      1 / exp(x$icoef["Log(scale)"])
    })
    eta0 <- lapply(seq_along(fit), function(g) {
      exp(fit[[g]]$coefficients)^(-kappa0[[g]])
    })

    # initials
    beta.ini <- rep(list(Beta.ini), S)
    gamma.ini <- rep(list(initial$gamma.ini), S)
    log.like <- lapply(survObj, function(x) {
      coxph(Surv(x$t, x$di, type = c("right")) ~ 1)$loglik
    })

    # survival object and covariate data
    survObj <- list(
      "t" = lapply(survObj, function(x) x$t),
      "di" = lapply(survObj, function(x) x$di),
      "X" = lapply(survObj, function(x) x$X), # each subgroup scaled separately
      "n" = lapply(survObj, function(x) length(x$t)) # sample size per subgroup
    )

    survObj$SSig <- lapply(survObj$X, function(x) t(x) %*% x) # sample covariance matrix per subgroup
    survObj$p <- ncol(survObj$X[[1]])
  }

  # check the formula
  cl <- match.call()

  # set hyperparamters of all piors

  if (model.type == "Sub-struct" && S > 1) {
    MRF2b <- TRUE
    # b02 = 0
    hyperpar$b <- hyperpar$b[1]
  }
  # if (MRF2b) {
  #   # b0 = c(b01, b02)
  #   b0 <- hyperpar$b
  # }

  # defining parameters for the main MCMC simulation function

  hyperpar$kappa0 <- kappa0
  hyperpar$eta0 <- eta0
  if (!MRF.G) {
    # hyperparameters for graphical learning
    hyperpar$pi.G <- 2 / (p - 1) # prior edge selection probability for all edges
    if (is.null(hyperpar$nu0)) {
      stop("Please specify hyperparameter 'hyperpar$nu0'!")
    }
    if (is.null(hyperpar$nu1)) {
      stop("Please specify hyperparameter 'hyperpar$nu1'!")
    }
    hyperpar$V0 <- hyperpar$nu0^2 * matrix(1, p, p) # matrix with small variances for prior of precision matrices
    hyperpar$V1 <- hyperpar$nu1^2 * matrix(1, p, p) # matrix with large variances for prior of precision matrices
  }
  if (!"be.prop.sd.scale" %in% names(hyperpar)) {
    if (S == 1 && MRF.G) {
      hyperpar$be.prop.sd.scale <- 2.4
    } else {
      hyperpar$be.prop.sd.scale <- rep(list(2.4), S)
    }
  }

  initial <- list("gamma.ini" = gamma.ini, "beta.ini" = beta.ini, "log.like.ini" = log.like)

  # if(model.type %in% c("CoxBVSSL", "Sub-struct") || (model.type == "Pooled" && !MRF.G)){
  if (!MRF.G) {
    # starting value for graph G (pS x pS matrix across all subgroups)
    if (g.ini %in% c(0, 1)) {
      G_ss <- matrix(g.ini, p, p) # subgraph within each subgroup
      diag(G_ss) <- 0

      G_rs <- diag(g.ini, p, p) # subgraph between subgroups
      if (model.type == "Sub-struct" || (model.type == "Pooled" && !MRF.G)) {
        G_rs <- matrix(0, p, p)
      }

      G.ini <- do.call("cbind", rep(list(do.call("rbind", rep(list(G_rs), S))), S)) # complete graph G
      for (g in 1:S) {
        G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <- G_ss
      }
      # variances for prior of precision matrices:
      if (g.ini == 0) {
        V.ini <- rep(list(hyperpar$V0), S)
      } else {
        V.ini <- rep(list(hyperpar$V1), S)
      }
    } else { # g.ini is rate of selected edges in subgraphs within each subgroup

      G_rs <- diag(1, p, p) # subgraph between subgroups (since all subgroups have the same starting values for gamma)
      G_ss <- matrix(0, p, p) # subgraph within each subgroup
      updiag <- sample(which(upper.tri(G_ss)), round(((p^2 - p) / 2) * g.ini))
      G_ss[updiag] <- 1
      G_ss[lower.tri(G_ss)] <- t(G_ss)[lower.tri(G_ss)]

      if (model.type == "Sub-struct" || (model.type == "Pooled" && !MRF.G)) {
        G_rs <- matrix(0, p, p)
      }

      G.ini <- do.call("cbind", rep(list(do.call("rbind", rep(list(G_rs), S))), S))
      for (g in 1:S) {
        G.ini[(g - 1) * p + (1:p), (g - 1) * p + (1:p)] <- G_ss
      }

      # variances for prior of precision matrices
      v_ini <- hyperpar$V1
      v_ini[which(G_ss == 0)] <- hyperpar$nu0^2
      V.ini <- rep(list(v_ini), S)
    }

    # list with precision and covariance matrix per subgroup
    C.ini <- Sig.ini <- rep(list(diag(1, p, p)), S)

    initial <- list(
      "G.ini" = G.ini, "V.ini" = V.ini, "C.ini" = C.ini, "Sig.ini" = Sig.ini,
      "gamma.ini" = gamma.ini, "beta.ini" = beta.ini, "log.like.ini" = log.like
    )
  }

  if (MRF.G) {
    initial$G.ini <- data.matrix(hyperpar$G)
    if (model.type == "Pooled" && any(dim(initial$G.ini) != c(p, p))) {
      stop("Hyperparameter 'hyperpar$G' has incorrect dimensions!")
    }
    if (model.type != "Pooled" && any(dim(initial$G.ini) != c(p * S, p * S))) {
      stop("Hyperparameter 'hyperpar$G' has incorrect dimensions!")
    }
  }

  ret <- list(input = list(), output = list(), call = cl)
  class(ret) <- "BayesSurvive"

  # Copy the inputs
  ret$input["p"] <- survObj$p
  ret$input["nIter"] <- nIter
  ret$input["burnin"] <- burnin
  ret$input["thin"] <- 1
  ret$input["rw"] <- NULL
  ret$input["S"] <- S
  ret$input["model.type"] <- model.type
  ret$input["MRF2b"] <- MRF2b
  ret$input["MRF.G"] <- MRF.G
  ret$input$hyperpar <- hyperpar

  ret$output <- func_MCMC(
    survObj = survObj,
    hyperpar = hyperpar,
    initial = initial,
    nIter = nIter,
    burnin = burnin,
    thin = thin,
    S = S,
    method = model.type,
    MRF_2b = MRF2b,
    MRF_G = MRF.G,
    output_graph_para,
    verbose,
    cpp
  )

  if (S == 1 && MRF.G) {
    ret$output$survObj <- list(
      "t" = survObj$t,
      "di" = survObj$di,
      "lp.mean" = NULL,
      "lp.all" = NULL
    )
    ret$output$survObj$lp.mean <- survObj$X %*%
      matrix(ret$output$beta.margin, ncol = 1)
    ret$output$survObj$lp.all <- survObj$X %*% t(ret$output$beta.p)
  } else {
    ret$output$survObj <- vector("list", S)
    for (s in 1:S) {
      ret$output$survObj[[s]]$t <- survObj$t[[s]]
      ret$output$survObj[[s]]$di <- survObj$di[[s]]
      ret$output$survObj[[s]]$lp.mean <-
        survObj$X[[s]] %*% matrix(ret$output$beta.margin[[1]], ncol = 1)

      ret$output$survObj[[s]]$lp.all <-
        survObj$X[[s]] %*% t(ret$output$beta.p[[s]])
    }
  }

  return(ret)
}
