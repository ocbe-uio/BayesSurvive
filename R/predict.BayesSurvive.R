#' @title Predict survival risk
#'
#' @description
#' Predict survival probability, (cumulative) hazard or (integrated) Brier
#' scores based on Cox regression models
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
#' @param verbose logical value to print IBS of the NULL model and the Bayesian
#' Cox model
#' @param \dots not used
#'
#' @return A list object with 5 components if \code{type="brier"} including
#' "model", "times", "Brier", "IBS" and "IPA" (Index of Prediction Accuracy),
#' otherwise  a list of 7 components with the first component as
#' the specified argument \code{type} and "se", "band", "type", "diag",
#' "baseline" and "times", see function \code{riskRegression::predictCox} for
#' details
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
#' predict(fit, survObj.new = dataset)
#' }
#'
#' @export
predict.BayesSurvive <- function(object, survObj.new, type = "brier",
                                 method = "mean", times = NULL, subgroup = 1,
                                 verbose = TRUE, ...) {
  if (!inherits(object, "BayesSurvive")) {
    stop("Use only with 'BayesSurvive' object!")
  }

  if (any(!type %in% c("brier", "hazard", "cumhazard", "survival")) && length(type) == 1) {
    stop("Argument 'type' has to be one of 'brier', 'hazard', 'cumhazard'' or 'survival'!")
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
    survObj$lp.all <- survObj$lp.all[, -1]
  } else {
    beta_m <- object$output$beta.margin[[subgroup]]
    survObj.new <- survObj.new[[subgroup]]
    betas <- object$output$beta.p[[subgroup]][-(1:(object$input$burnin / object$input$thin + 1)), ]
    survObj <- object$output$survObj[[subgroup]]
    survObj$lp.all <- survObj$lp.all[, -(1:(object$input$burnin / object$input$thin + 1))]
  }

  ibs <- data.frame("Null model" = rep(NA, 3), "Bayesian Cox" = rep(NA, 3))

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
    return(fit)
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
          metrics = "brier", summary = c("ibs", "ipa"),
          null.model = TRUE, times = times
        )$Brier$score

        # Brier <- BrierScore[BrierScore$model != "Null model", ]
        # extract scores for Null model and the Bayesian Cox model
        ibs[1, ] <- Brier$Brier[c(length(times), length(times) * 2)]
        ibs[2, ] <- Brier$IBS[c(length(times), length(times) * 2)]
        ibs[3, ] <- Brier$IPA[c(length(times), length(times) * 2)]
      }
    } else {
      # obtain linear predictors corresponding to MCMC estimates
      # lp_all_train <- survObj_train[, -c(1:2)] %*% t(betas)
      lp_all_test <- survObj.new$X %*% t(betas)
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
      BrierScore <- Score(list("Bayesian Cox" = model_train),
        formula = Surv(time, status) ~ 1,
        data = data_test, conf.int = FALSE,
        metrics = "brier", summary = c("ibs", "ipa"),
        null.model = TRUE, times = times
      )$Brier$score
      BrierScore[is.na(BrierScore)] <- 0
      # Brier scores of NULL.model do not change
      Brier.null <- BrierScore[1:(nrow(BrierScore) / 2), -c(1:2)]
      Brier <- BrierScore[-c(1:(nrow(BrierScore) / 2)), -c(1:2)]

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
          metrics = "brier", summary = c("ibs", "ipa"),
          null.model = TRUE, times = times
        )$Brier$score

        BrierScore[is.na(BrierScore)] <- 0
        Brier <- Brier + BrierScore[-c(1:(nrow(BrierScore) / 2)), -c(1:2)]
      }
      Brier <- rbind(Brier.null, Brier / nrow(betas))

      # extract IBS for Null model and the Bayesian Cox model
      ibs[1, ] <- Brier$Brier[c(length(times), length(times) * 2)]
      ibs[2, ] <- Brier$IBS[c(length(times), length(times) * 2)]
      ibs[3, ] <- Brier$IPA[c(length(times), length(times) * 2)]
    }

    # ibs <- rbind(rep(max(times), 2), ibs)
    t0 <- round(max(times), digits = 1)
    rownames(ibs) <- c(
      paste0("Brier(t=", t0, ")"),
      paste0("IBS(t:0~", t0, ")"),
      paste0("IPA(t=", t0, ")")
    )
    if (verbose) {
      # cat("                      IBS\n",
      #   "  Null model          ", ibs[1, 1],
      #   "\n  Bayesian Cox model  ", ibs[1, 2], "\n",
      #   sep = ""
      # )
      print(t(ibs))
    }
    invisible(Brier)
    # return(Brier)
  }
}
