#' @title Predict survival risk
#'
#' @description
#' Predict survival probability, (cumulative) hazard or (integrated) Brier scores based on Cox regression models
#'
#' @name predict.BayesSurvive
#'
#' @importFrom riskRegression Score predictCox
#' @importFrom survival coxph Surv
#'
#' @param object fitted object obtained with \code{BayesSurvive}
#' @param survObj.new a list containing observed data from new subjects with
#' components \code{t}, \code{di}, \code{x}. If \code{type} is among
#' \code{c("hazard", "cumhazard", "survival")}, only \code{survObj.new$X} is
#' needed.
#' @param times time points at which to evaluate the risks. If \code{NULL}
#' (default), the event/censoring times are used. If \code{type="brier"}, the
#' largest one of the \code{times} is used
#' @param type option to chose for predicting brier scores with \code{type="brier"}
#' or one of \code{type=c("brier", "hazard", "cumhazard", "survival")})
#' @param method option to use the posterior mean (\code{"mean"}) of coefficients
#' for prediction or Bayesian model averaging (\code{"BMA"}) for prediction
#' @param subgroup index of the subgroup in \code{survObj.new} for prediction.
#' Default value is 1
#' @param \dots not used
#'
#' @keywords survival
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
#' predict(fit, survObj.new = dataset)
#' }
#'
#' @export
predict.BayesSurvive <- function(object, survObj.new, type = "brier",
                              method = "mean", times = NULL, subgroup = 1,
                              ...) {
  if (!inherits(object, "BayesSurvive")) {
    stop("Use only with 'BayesSurvive' object!")
  }

  if (any(!type %in% c("brier", "hazard", "cumhazard", "survival"))) {
    stop("Argument 'type' has to be among 'brier', 'hazard', 'cumhazard'' or 'survival'!")
  }

  if (!method %in% c("mean", "BMA")) {
    stop("Argument 'method' has to be 'mean' or 'BMA'!")
  }

  # Check the survival input object
  if (is.null(survObj.new)) {
    stop("Data 'survObj.new' must be provided!")
  }

  # survObj_train <- as.matrix(read.table(paste0(object$input$outFilePath, "data.txt")))
  if (object$input$S == 1 && object$input$MRF.G) {
    beta_m <- object$output$beta.margin # colMeans(betas)
    betas <- object$output$beta.p[-(1:(object$input$burnin / object$input$thin + 1)), ]
    survObj <- object$output$survObj
  } else {
    beta_m <- object$output$beta.margin[[subgroup]]
    survObj.new <- survObj.new[[subgroup]]
    betas <- object$output$beta.p[[subgroup]][-(1:(object$input$burnin / object$input$thin + 1)), ]
    survObj <- object$output$survObj[[subgroup]]
    survObj$lp.all <- survObj$lp.all[, -(1:(object$input$burnin / object$input$thin + 1))]
  }

  ibs <- data.frame("Null model" = NA, "Bayesian Cox" = NA)

  if (is.null(times)) {
    times <- sort(unique(survObj.new$t))
  }
  if (max(times) > max(survObj$t)) {
    message("NOTE: The evaluation times were truncated by the largest observation time in the training data!")
    times <- times[times <= max(survObj$t)]
  }

  if (length(type) == 1 && type == "brier") {
    times <- seq(0, max(times), length = 100)
  }

  if (!"brier" %in% type) {
    if (method == "mean") {
      # lp <- as.vector(survObj[, -c(1:2)] %*% beta_m)
      data_train <- data.frame(
        time = survObj$t,
        status = survObj$di, lp = survObj$lp.mean
      )
      model_train <- coxph(Surv(time, status) ~ lp,
        data = data_train, y = TRUE, x = TRUE
      )

      lp_test <- as.vector(survObj.new$X %*% beta_m)
      # data_test <- data.frame(time = survObj.new$t, status = survObj.new$di,
      #                         lp = lp_test)
      fit <- predictCox(model_train,
        times = times, type = type,
        newdata = data.frame(lp = lp_test)
      )
      fit[names(fit) %in% c(
        "lastEventTime", "nTimes",
        "var.lp", "var.strata"
      )] <- NULL
    } else {
      # obtain linear predictors corresponding to MCMC estimates
      # lp_all_train <- survObj[, -c(1:2)] %*% t(betas)
      lp_all_test <- survObj.new$X %*% t(betas)
      data_train <- data.frame(
        time = survObj$t,
        status = survObj$di,
        lp = survObj$lp.mean[, 1]
      )
      # data_test <- data.frame(time = survObj.new$t,
      #                         status = survObj.new$di,
      #                         lp = lp_all_test[, 1])
      model_train <- coxph(Surv(time, status) ~ lp,
        data = data_train,
        y = TRUE, x = TRUE
      )
      # calculate survival prediction based on the 1st MCMC estimates
      fit <- predictCox(model_train,
        times = times, type = type,
        newdata = data.frame(lp = lp_all_test[, 1])
      )
      fit0 <- fit[names(fit) %in% type]

      # calculate Brier scores based on other MCMC estimates
      for (i in 2:nrow(betas)) {
        data_train$lp <- survObj$lp.mean[, i]
        # data_test$lp <- lp_all_test[, i]
        model_train <- coxph(Surv(time, status) ~ lp,
          data = data_train,
          y = TRUE, x = TRUE
        )
        fit <- predictCox(model_train,
          times = times, type = type,
          newdata = data.frame(lp = lp_all_test[, i])
        )
        fit0 <- Map("+", fit0, fit[names(fit) %in% type])
      }
      fit0 <- Map("/", fit0, nrow(betas))
      fit[names(fit) %in% type] <- fit0

      fit[names(fit) %in% c(
        "lastEventTime", "nTimes",
        "var.lp", "var.strata"
      )] <- NULL
    }
    fit
  } else {
    if (method == "mean") {
      # beta_m <- colMeans(betas)
      # lp <- as.vector(survObj_train[, -c(1:2)] %*% beta_m)
      data_train <- data.frame(
        time = survObj$t,
        status = survObj$di,
        lp = survObj$lp.mean
      )
      model_train <- coxph(Surv(time, status) ~ lp, data = data_train, y = TRUE, x = TRUE)

      lp_test <- as.vector(survObj.new$X %*% beta_m)
      data_test <- data.frame(time = survObj.new$t, status = survObj.new$di, lp = lp_test)

      if (any(type %in% c("hazard", "cumhazard", "survival"))) {
        Brier <- predictCox(model_train,
          times = times, type = type,
          newdata = data.frame(survObj.new$X)
        )
        Brier <- data.frame(times = fit$times, fit[names(fit) %in% type])
      } else {
        Brier <- Score(list("Bayesian Cox" = model_train),
          formula = Surv(time, status) ~ 1,
          data = data_test, conf.int = FALSE,
          metrics = "brier", summary = "ibs",
          times = times
        )$Brier$score
        # Brier <- BrierScore[BrierScore$model != "Null model", ]
        # extract IBS for Null model and the Bayesian Cox model
        ibs[1, ] <- Brier$IBS[c(length(times), length(times) * 2)]
      }
    } else {
      # obtain linear predictors corresponding to MCMC estimates
      # lp_all_train <- survObj_train[, -c(1:2)] %*% t(betas)
      lp_all_test <- survObj.new$X %*% t(betas)
      Brier <- matrix(0, nrow = length(times) * 2, ncol = 2)
      data_train <- data.frame(
        time = survObj$t,
        status = survObj$di,
        lp = survObj$lp.all[, 1]
      )
      data_test <- data.frame(
        time = survObj.new$t,
        status = survObj.new$di,
        lp = lp_all_test[, 1]
      )
      model_train <- coxph(Surv(time, status) ~ lp,
        data = data_train,
        y = TRUE, x = TRUE
      )
      # calculate Brier scores based on the 1st MCMC estimates
      BrierScore0 <- Score(list("Bayesian Cox" = model_train),
        formula = Surv(time, status) ~ 1,
        data = data_test, conf.int = FALSE,
        metrics = "brier", summary = "ibs",
        times = times
      )$Brier$score
      Brier <- as.matrix(BrierScore0[1:(nrow(BrierScore0) / 2), -c(1:2)])
      # calculate Brier scores based on other MCMC estimates
      for (i in 2:nrow(betas)) {
        data_train$lp <- survObj$lp.all[, i] # lp_all_train[, i]
        data_test$lp <- lp_all_test[, i]
        model_train <- coxph(Surv(time, status) ~ lp,
          data = data_train,
          y = TRUE, x = TRUE
        )
        BrierScore <- Score(list("Bayesian Cox" = model_train),
          formula = Surv(time, status) ~ 1,
          data = data_test, conf.int = FALSE,
          metrics = "brier", summary = "ibs",
          null.model = FALSE, times = times
        )$Brier$score
        Brier <- Brier + as.matrix(BrierScore[, -c(1:2)])
      }
      BrierScore0 <- data.frame(BrierScore0)
      BrierScore0[-c(1:(nrow(BrierScore0) / 2)), 3:4] <- Brier / nrow(betas)
      Brier <- BrierScore0
      # extract IBS for Null model and the Bayesian Cox model
      ibs[1, ] <- Brier$IBS[c(length(times), length(times) * 2)]
    }
    rownames(ibs) <- "IBS"
    print(ibs)
    invisible(Brier)
  }

  # invisible(fit)
}
