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

  if (x$input$S > 1 || !x$input$MRF.G) {
    x$output$beta.p <- x$output$beta.p[[subgroup]]
  }

  if (inherits(x, "BayesSurvive")) {
    if (is.null(colnames(x$output$beta.p))) {
      x_names <- paste0("x", seq_len(ncol(x$output$beta.p)))
    } else {
      x_names <- colnames(x$output$beta.p)
    }
    beta_p <- x$output$beta.p[-(1:(x$input$burnin / x$input$thin + 1)), ]
  } else {
    if (is.null(colnames(x))) {
      x_names <- paste0("x", seq_len(ncol(x)))
    } else {
      x_names <- colnames(x)
    }
    beta_p <- x
  }

  # pdf("psbcBeta.pdf", height = 5, width = 3.5)
  beta_est <- apply(beta_p, 2, type)
  beta_L <- apply(beta_p, 2, quantile, 0.025)
  beta_U <- apply(beta_p, 2, quantile, 0.975)
  tbl <- data.frame(term = x_names, estimate = beta_est, conf.low = beta_L, conf.high = beta_U)
  tbl$term <- factor(tbl$term, levels = tbl$term)

  # Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")
  pCoef <- ggcoef(tbl, ...) + xlab(expression(Posterior ~ ~beta)) + ylab("")
  pCoef
  # dev.off()
}
