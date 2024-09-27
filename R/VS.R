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
#' @param x fitted object obtained with \code{BayesSurvive}, or a matrix/array,
#' or a list consisting of matrices and arrays
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
  if (!(inherits(x, "BayesSurvive") || is.array(x) || is.list(x))) {
    stop("Use only with 'BayesSurvive' object or a matrix or an array or a list!")
  }
  if (!inherits(x, "BayesSurvive") && is.list(x)) {
    if (any(!sapply(x, function(xx) is.array(xx)))) {
      stop("The list input has to consist of matrices and/or arrays!")
    }
    # If the input is a matrix or array, reformat it to be an list
    if (!is.list(x)) {
      x <- list(x)
    }
  }
  if (!method %in% c("CI", "SNC", "MPM", "FDR")) { # "SNC-BIC",
    stop("'method' should be one of c('CI', 'SNC', 'MPM', 'FDR')!")
  }

  if (method == "CI") {
    # for an input from "BayesSurvive"
    if (inherits(x, "BayesSurvive")) {
      betas <- coef.BayesSurvive(x,
        type = "mean", CI = 95,
        subgroup = subgroup
      )
      ret <- (betas$CI.lower > 0) | (betas$CI.upper < 0)
    } else {
      # for an output from a matrix or array or list
      if (!is.list(x)) {
        betas <- list()
        # if (is.matrix(x)) {
        #   betas$CI.lower <- apply(x, 2, quantile, 0.025)
        #   betas$CI.upper <- apply(x, 2, quantile, 0.975)
        # }
        # if (is.array(x)) {
        #   betas$CI.lower <- apply(x, c(2,3), quantile, 0.025)
        #   betas$CI.upper <- apply(x, c(2,3), quantile, 0.975)
        # }
        # ret <- (betas$CI.lower > 0) | (betas$CI.upper < 0)
      } else {
        betas <- rep(list(NULL), length(x))
        ret <- rep(list(NULL), length(x))
        for (l in seq_len(length(x))) {
          # if (length(dim(x[[l]])) == 2) {
          #   betas[[l]]$CI.lower <- apply(x[[l]], 2, quantile, 0.025)
          #   betas[[l]]$CI.upper <- apply(x[[l]], 2, quantile, 0.975)
          # }
          # if (length(dim(x[[l]])) == 3) {
          #   betas[[l]]$CI.lower <- apply(x[[l]], c(2,3), quantile, 0.025)
          #   betas[[l]]$CI.upper <- apply(x[[l]], c(2,3), quantile, 0.975)
          # }
          betas[[l]]$CI.lower <- apply(x[[l]], seq_len(length(dim(x[[l]])))[-1], quantile, 0.025)
          betas[[l]]$CI.upper <- apply(x[[l]], seq_len(length(dim(x[[l]])))[-1], quantile, 0.975)
          ret[[l]] <- (betas[[l]]$CI.lower > 0) | (betas[[l]]$CI.upper < 0)
        }
      }
    }
  }

  if (method == "SNC") {
    if (is.na(threshold)) {
      threshold <- 0.5
    }

    if (inherits(x, "BayesSurvive")) {
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
    } else {
      # count the total number of parameters in the list
      total_num <- 0
      for (l in seq_len(length(x))) {
        total_num <- total_num + prod(dim(x[[l]])[-1])
      }

      ret <- rep(list(NULL), length(x))
      for (l in seq_len(length(x))) {
        # define the size of each component in the list
        ret[[l]] <- array(FALSE, dim = dim(x[[l]])[-1])

        if (is.matrix(x[[l]])) { # for an matrix
          for (j in seq_len(NCOL(x[[l]]))) {
            if (sum(abs(x[[l]][, j]) > sd(x[[l]][, j])) / total_num > threshold) {
              ret[[l]][j] <- TRUE
            }
          }
        } else { # for an array
          for (j in seq_len(dim(x[[l]])[2])) {
            for (k in seq_len(dim(x[[l]])[3])) {
              if (sum(abs(x[[l]][, j, k]) > sd(x[[l]][, j, k])) / total_num > threshold) {
                ret[[l]][j, k] <- TRUE
              }
            }
          }
        }
      }
    }
  }

  if (method == "MPM") {
    if (inherits(x, "BayesSurvive")) {
      if (x$input$S > 1 || !x$input$MRF.G) {
        x$output$gamma.margin <- x$output$gamma.margin[[subgroup]]
      }
      ret <- (x$output$gamma.margin >= 0.5)
    } else {
      ret <- rep(list(NULL), length(x))

      for (l in seq_len(length(x))) {
        # define the size of each component in the list
        ret[[l]] <- array(FALSE, dim = dim(x[[l]])[-1])

        if (is.matrix(x[[l]])) { # for an matrix
          for (j in seq_len(NCOL(x[[l]]))) {
            ret[[l]][j] <- (mean(x[[l]][, j] != 0) >= 0.5)
          }
        } else { # for an array
          for (j in seq_len(dim(x[[l]])[2])) {
            for (k in seq_len(dim(x[[l]])[3])) {
              ret[[l]][j, k] <- (mean(x[[l]][, j, k] != 0) >= 0.5)
            }
          }
        }
      }
    }
  }

  if (method == "FDR") {
    if (is.na(threshold)) {
      threshold <- 0.05
    }
    if (inherits(x, "BayesSurvive")) {
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
    } else {
      ret <- gammas <- rep(list(NULL), length(x))
      gammas_vec <- NULL

      # compute mPIPs
      for (l in seq_len(length(x))) {
        # define the size of each component in the list
        gammas[[l]] <- array(FALSE, dim = dim(x[[l]])[-1])

        if (is.matrix(x[[l]])) { # for an matrix
          for (j in seq_len(NCOL(x[[l]]))) {
            gammas[[l]][j] <- mean(x[[l]][, j] != 0)
          }
        } else { # for an array
          for (j in seq_len(dim(x[[l]])[2])) {
            for (k in seq_len(dim(x[[l]])[3])) {
              gammas[[l]][j, k] <- mean(x[[l]][, j, k] != 0)
            }
          }
        }
        # save all mPIPs into a vector
        gammas_vec <- c(gammas_vec, as.vector(gammas[[l]]))
      }

      sorted_gammas <- sort(gammas_vec, decreasing = TRUE)
      # computing the fdr
      fdr <- cumsum((1 - sorted_gammas)) / seq_len(length(sorted_gammas))
      # determine index of the largest fdr less than threshold
      if (min(fdr) >= threshold) {
        ret <- rep(FALSE, length(gammas_vec))
      } else {
        thecut.index <- max(which(fdr < threshold))
        gammas_threshold <- sorted_gammas[thecut.index]
        ret_vec <- (gammas_vec > gammas_threshold)
      }

      # reformat the results into a list
      for (l in seq_len(length(x))) {
        ret[[l]] <- array(FALSE, dim = dim(x[[l]])[-1])
        len <- prod(dim(x[[l]])[-1])
        ret[[l]] <- array(ret_vec[1:len], dim = dim(x[[l]])[-1])
        ret_vec <- ret_vec[-c(1:len)]
      }
    }
  }

  ret
}
