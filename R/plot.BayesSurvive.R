#' @title Create a plot of estimated coefficients
#' @description
#' Plot point estimates of regression coefficients and 95\% credible intervals
#'
#' @name plot.BayesSurvive
#'
#' @importFrom GGally ggcoef
#' @importFrom ggplot2 xlab ylab
#' @importFrom stats quantile
#'
#' @param x an object of class \code{BayesSurvive} or a matrix. If \code{x}
#' is a matrix, use \code{BayesSurvive:::plot.BayesSurvive(x)}
#' @param type type of point estimates of regression coefficients. One of
#' \code{c("mean", "median")}. Default is \code{mean}
#' @param interval logical argument to show 95\% credible intervals. Default
#' is \code{TRUE}
#' @param subgroup index of the subgroup for visualizing posterior coefficients
#' @param ... additional arguments sent to \code{ggplot2::geom_point()}
#'
#' @return ggplot object
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
#' # show posterior mean of coefficients and 95% credible intervals
#' library("GGally")
#' plot(fit) +
#'   coord_flip() +
#'   theme(axis.text.x = element_text(angle = 90, size = 7))
#' }
#'
#' @export
plot.BayesSurvive <- function(x, type = "mean", interval = TRUE,
                              subgroup = 1, ...) {
  if (!(inherits(x, "BayesSurvive") || is.matrix(x))) {
    stop("Use only with 'BayesSurvive' object or a matrix!")
  }

  if (length(type) == 1) {
    if (!type %in% c("mean", "median")) {
      stop("'type' should be one of c('mean', 'median')!")
    }
  } else {
    stop("'type' should be one of c('mean', 'median')!")
  }

  if (!is.logical(interval)) {
    stop("Argument 'interval' must be a logical value!")
  }

  if (inherits(x, "BayesSurvive")) {
    tbl <- coef.BayesSurvive(x,
      type = type, CI = 95,
      subgroup = subgroup
    )
    names(tbl)[2:4] <- c("estimate", "conf.low", "conf.high")
    tbl$term <- factor(tbl$term, levels = tbl$term)
  } else {
    if (is.null(colnames(x))) {
      x_names <- paste0("x", seq_len(ncol(x)))
    } else {
      x_names <- colnames(x)
    }
    beta_p <- x

    beta_est <- apply(beta_p, 2, type)
    beta_L <- apply(beta_p, 2, quantile, 0.025)
    beta_U <- apply(beta_p, 2, quantile, 0.975)
    tbl <- data.frame(
      term = x_names, estimate = beta_est,
      conf.low = beta_L, conf.high = beta_U
    )
    tbl$term <- factor(tbl$term, levels = tbl$term)
  }

  # pdf("psbcBeta.pdf", height = 5, width = 3.5)

  # Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")
  pCoef <- ggcoef(tbl, conf.int = interval, ...) +
    xlab(expression(Posterior ~ ~beta)) + ylab("")
  pCoef
  # dev.off()
}
