#' @title Time-dependent Brier scores
#'
#' @description
#' Predict time-dependent Brier scores based on Cox regression models
#'
#' @name plotBrier
#'
#' @importFrom ggplot2 ggplot aes geom_step theme element_blank
#'
#' @param object fitted object obtained with \code{BayesSurvive}
#' @param survObj.new a list containing observed data from new subjects with
#' components \code{t}, \code{di}, \code{X}
#' @param times maximum time point to evaluate the prediction
#' @param method option to use the posterior mean (\code{"mean"}) of coefficients
#' for prediction or Bayesian model averaging (\code{"BMA"}) for prediction
#' @param subgroup index of the subgroup in \code{survObj.new} for prediction.
#' Default value is 1
#' 
#' @return A \code{ggplot2::ggplot} object. See \code{?ggplot2::ggplot} for more
#' details of the object. 
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
#' # predict survival probabilities of the train data
#' plotBrier(fit, survObj.new = dataset)
#' }
#'
#' @export
plotBrier <- function(object, survObj.new = NULL,
                      method = "mean", times = NULL, subgroup = 1) {
  # if (is.null(times))
  #   times <- sort(unique(survObj.new$t))
  capture.output(
    Brier_score <- predict.BayesSurvive(object,
      survObj.new = survObj.new,
      method = method,
      times = times,
      subgroup = subgroup
    ),
    file = nullfile()
  )

  Brier <- model <- NULL
  # Brier_score %>%
  ggplot2::ggplot(
    Brier_score,
    aes(times, Brier, group = model, color = model)
  ) +
    xlab("Evaluation time points") +
    ylab("Brier score") +
    geom_step(direction = "vh") +
    theme(legend.title = element_blank())
}
