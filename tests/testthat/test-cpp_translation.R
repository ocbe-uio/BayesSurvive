# Load the example data
data <- list(
  "X" = simData[[1]]$X,
  "t" = simData[[1]]$time,
  "di" = simData[[1]]$status
)
data_2S <- list(
  data,
  list(
    "X" = simData[[2]]$X,
    "t" = simData[[2]]$time,
    "di" = simData[[2]]$status
  )
)
# Run a Bayesian Cox model

## Initial value: null model without covariates
initial <- list("gamma.ini" = rep(0, ncol(data$X)))

# Prior parameters
hyperPooled = list(
  "c0"     = 2,                      # prior of baseline hazard
  "tau"    = 0.0375,                 # sd (spike) for coefficient prior
  "cb"     = 20,                     # sd (slab) for coefficient prior
  "pi.ga"  = 0.02,                   # prior variable selection probability for standard Cox models
  "a"      = -4,                     # hyperparameter in MRF prior
  "b"      = 0.1,                    # hyperparameter in MRF prior
  "G"      = simData$G               # hyperparameter in MRF prior
)
hyperPooled_2S <- hyperPooled
hyperPooled_2S$G <- Matrix::bdiag(simData$G, simData$G)

# Run a 'Pooled' Bayesian Cox model with graphical learning

set.seed(715074)
BayesSurvive_wrap <- function(
  data, initial, hyper, model = "Pooled", use_cpp = FALSE, n_iter = 5,
  MRF_G = TRUE, MRF_2b = FALSE, verbose = FALSE
  ) {
  if (!MRF_G) {
    if (!is.null(names(data))) {
      data <- list(data)
      # hyper$b <- c(0.1, 0.2) # TODO: uncomment for !MRF_2b cases
    }
    hyper$lambda <- 3 # TODO: mandatory for !MRG.G? Add validation!
    hyper$nu0 <- 0.05
    hyper$nu1 <- 5
  }
  BayesSurvive(
    survObj = data, model.type = model, MRF.G = MRF_G, MRF2b = MRF_2b,
    verbose = verbose, hyperpar = hyper, initial = initial, nIter = n_iter,
    burnin = floor(n_iter / 2), cpp = use_cpp
  )
}
fit_R <- BayesSurvive_wrap(data, initial, hyperPooled)
fit_C <- BayesSurvive_wrap(data, initial, hyperPooled, use_cpp = TRUE)
fit_R2S <- BayesSurvive_wrap(data_2S, initial, hyperPooled_2S, "CoxBVSSL")
fit_C2S <- BayesSurvive_wrap(data_2S, initial, hyperPooled_2S, "CoxBVSSL", use_cpp = TRUE)
fit_R_noMRFG <- BayesSurvive_wrap(data, initial, hyperPooled, MRF_G = FALSE, n_iter = 2L)
fit_C_noMRFG <- BayesSurvive_wrap(data, initial, hyperPooled, MRF_G = FALSE, use_cpp = TRUE, n_iter = 2L) # FIXME: Error: dot(): objects must have the same number of elements
fit_R_2b <- BayesSurvive_wrap(data, initial, hyperPooled, MRF_2b = TRUE)
fit_C_2b <- BayesSurvive_wrap(data, initial, hyperPooled, MRF_2b = TRUE, use_cpp = TRUE)
fit_R_2b_no_G <- BayesSurvive_wrap(data_2S, initial, hyperPooled_2S, MRF_2b = TRUE, MRF_G = FALSE, n_iter = 2L)
fit_C_2b_no_G <- BayesSurvive_wrap(data_2S, initial, hyperPooled_2S, MRF_2b = TRUE, MRF_G = FALSE, use_cpp = TRUE, n_iter = 2L) # FIXME: Error: Mat::operator(): index out of bounds

# TODO: reduce. Takes 4 minutes!
# TODO: reorganize tests so that they come right after each fit_R/fit_C pair
test_that("R and C++ objects are similar", {
  expect_equal(fit_R$call, fit_C$call)
  expect_equal(fit_R$input, fit_C$input)
  for (obj in names(fit_R$output)[2]) {
    expect_equal(fit_R$output[[obj]], fit_C$output[[obj]], tolerance = 1)
  }
  expect_equal(fit_R2S$call, fit_C2S$call)
  expect_equal(fit_R2S$input, fit_C2S$input)
  for (obj in names(fit_R2S$output)[2]) {
    expect_equal(fit_R2S$output[[obj]], fit_C2S$output[[obj]], tolerance = 1)
  }
  expect_equal(fit_R_noMRFG$call, fit_C_noMRFG$call)
  expect_equal(fit_R_noMRFG$input, fit_C_noMRFG$input)
  for (obj in names(fit_R_noMRFG$output)[2]) {
    expect_equal(fit_R_noMRFG$output[[obj]], fit_C_noMRFG$output[[obj]], tolerance = 1)
  }
})
