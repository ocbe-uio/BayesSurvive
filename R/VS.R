#' @title Function to perform variable selection
#'
#' @description
#' Perform variable selection using the 95% credible interval (CI), scaled
#' neighborhood criterion (SNC), median probability model (MPM) or Bayesian
#' false discovery rate (FDR). Note that the Bayesian FDR only applies for each
#' subgroup if there are subgroups.
#'
#' @name VS
#'
#' @param x fitted object obtained with \code{BayesSurvive}
#' @param method variable selection method to choose from
#' \code{c("CI", "SNC", "MPM", "FDR")}. Default is "FDR"
#' @param threshold SNC threshold value (default 0.5) or the Bayesian expected false
#' discovery rate threshold (default 0.05)
#' @param subgroup index of the subgroup for visualizing posterior coefficients
#'
#' @return A boolean vector of selected (= TRUE) and rejected (= FALSE)
#' variables
#'
#' @references
#' Lee KH, Chakraborty S, Sun J (2015). Survival prediction and variable
#' selection with simultaneous shrinkage and grouping priors. Statistical
#' Analysis and Data Mining, 8:114-127
#'
#' Newton MA, Noueiry A, Sarkar D, Ahlquist P (2004). Detecting differential
#' gene expression with a semiparametric hierarchical mixture method.
#' Biostatistics, 5(2), 155-76
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
#' # show variable selection
#' VS(fit, method = "FDR")
#' }
#'
#' @export
VS <- function(x, method = "FDR", threshold = NA, subgroup = 1) {
  if (!inherits(x, "BayesSurvive")) {
    stop("Use only with 'BayesSurvive' object!")
  }
  if (!method %in% c("CI", "SNC", "SNC-BIC", "MPM", "FDR")) {
    stop("'method' should be one of c('CI', 'SNC', 'SNC-BIC', 'MPM', 'FDR')!")
  }

  if (method == "CI") {
    betas <- coef.BayesSurvive(x,
      type = "mean", CI = 95,
      subgroup = subgroup
    )
    ret <- (betas$CI.lower > 0) | (betas$CI.upper < 0)
  }

  if (method == "SNC") {
    if (is.na(threshold)) {
      threshold <- 0.5
    }
    if (x$input$S > 1 || !x$input$MRF.G) {
      x$output$beta.p <- x$output$beta.p[[subgroup]]
    } 
    beta_p <- x$output$beta.p[-(1:(x$input$burnin / x$input$thin + 1)), ]

    ret <- rep(FALSE, NCOL(beta_p))

    for (j in seq_len(NCOL(beta_p))) {
      if (sum(abs(beta_p[, j]) > sd(beta_p[, j])) / NROW(beta_p) > threshold) {
        ret[j] <- TRUE
      }
    }
  }

  if (method == "MPM") {
    if (x$input$S > 1 || !x$input$MRF.G) {
      x$output$gamma.margin <- x$output$gamma.margin[[subgroup]]
    }
    ret <- (x$output$gamma.margin >= 0.5)
  }

  if (method == "FDR") {
    if (is.na(threshold)) {
      threshold <- 0.05
    }
    if (x$input$S > 1 || !x$input$MRF.G) {
      x$output$gamma.margin <- x$output$gamma.margin[[subgroup]]
    }
    gammas <- x$output$gamma.margin

    sorted_gammas <- sort(gammas, decreasing = TRUE)
    # computing the fdr
    fdr <- cumsum((1 - sorted_gammas)) / seq_len(length(sorted_gammas))
    # determine index of the largest fdr less than threshold
    if (min(fdr) >= threshold) {
      ret <- rep(FALSE, length(gammas))
    } else {
      thecut.index <- max(which(fdr < threshold))
      gammas_threshold <- sorted_gammas[thecut.index]
      ret <- (gammas > gammas_threshold)
    }
  }

  ret
}
