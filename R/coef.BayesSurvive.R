#' @title Create a dataframe of estimated coefficients
#'
#' @description
#' Estimate regression coefficients with posterior mean/median, credible
#' intervals, standard deviation, or MPM estimates, posterior gammas
#'
#'
#' @name coef.BayesSurvive
#'
#' @importFrom stats quantile
#'
#' @param object an object of class \code{BayesSurvive}
#' @param MPM logical value to obtain MPM coefficients. Default: FALSE
#' @param type type of point estimates of regression coefficients. One of
#' \code{c("mean", "median")}. Default is \code{mean}
#' @param CI size (level, as a percentage) of the credible interval to report.
#' Default: 95, i.e. a 95\% credible interval
#' @param SD logical value to show each coefficient's standard deviation over
#' MCMC iterations
#' @param subgroup index of the subgroup for visualizing posterior coefficients
#' @param ... other arguments
#'
#' @return dataframe object
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
#'   "tau"    = 0.0375, # sd for coefficient prior
#'   "cb"     = 20, # sd for coefficient prior
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
#' # show posterior coefficients
#' betas <- coef(fit)
#' head(betas)
#' }
#'
#' @export
coef.BayesSurvive <- function(object, MPM = FALSE, type = "mean", CI = 95,
                              SD = FALSE, subgroup = 1, ...) {
  if (!(inherits(object, "BayesSurvive"))) {
    stop("Use only with 'BayesSurvive' object!")
  }

  if (length(type) == 1) {
    if (!type %in% c("mean", "median")) {
      stop("'type' should be one of c('mean', 'median')!")
    }
  } else {
    stop("'type' should be one of c('mean', 'median')!")
  }

  if (!(is.numeric(CI) && CI > 0 && CI < 100)) {
    stop("Argument 'CI' must be a numeric value between 0 and 100!")
  }
  if (!is.logical(SD)) {
    stop("Argument 'SD' must be logical!")
  }

  if (!MPM) {
    # posterior mean/median of coefficients
    if (object$input$S > 1 || !object$input$MRF.G) {
      object$output$beta.p <- object$output$beta.p[[subgroup]]
    }
    if (is.null(colnames(object$output$beta.p))) {
      x_names <- paste0("object", seq_len(ncol(object$output$beta.p)))
    } else {
      x_names <- colnames(object$output$beta.p)
    }
    beta_p <- object$output$beta.p[-(1:(object$input$burnin / object$input$thin + 1)), ]

    # pdf("psbcBeta.pdf", height = 5, width = 3.5)
    beta_est <- apply(beta_p, 2, type)
    beta_L <- apply(beta_p, 2, quantile, (1 - 0.01 * CI) / 2)
    beta_U <- apply(beta_p, 2, quantile, 0.5 + 0.01 * CI / 2)
    tbl <- data.frame(
      term = x_names, estimate = beta_est,
      CI.lower = beta_L, CI.upper = beta_U
    )
    names(tbl)[2] <- type
    if (SD) tbl <- data.frame(tbl, SD = apply(beta_p, 2, sd))
    tbl$term <- factor(tbl$term, levels = tbl$term)
  } else {
    # MPM coefficients
    if (object$input$S > 1 || !object$input$MRF.G) {
      object$output$beta.margin <- object$output$beta.margin[[subgroup]]
      object$output$gamma.margin <- object$output$gamma.margin[[subgroup]]
    }
    if (is.null(names(object$output$beta.margin))) {
      x_names <- paste0("object", seq_len(length(object$output$beta.margin)))
    } else {
      x_names <- names(object$output$beta.margin)
    }
    tbl <- data.frame(
      term = x_names,
      beta_MPM = object$output$beta.margin,
      gamma = object$output$gamma.margin
    )
    tbl$estimate <- (object$output$gamma.margin >= 0.5) *
      object$output$beta.margin / object$output$gamma.margin
    tbl$estimate[is.na(tbl$estimate)] <- 0
  }

  tbl
  
}
