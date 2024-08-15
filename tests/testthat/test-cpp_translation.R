# Load the example dataset
dataset <- list(
  "X" = simData[[1]]$X,
  "t" = simData[[1]]$time,
  "di" = simData[[1]]$status
)

# Run a Bayesian Cox model

## Initial value: null model without covariates
initial <- list("gamma.ini" = rep(0, ncol(dataset$X)))

# Prior parameters
hyperparPooled = list(
  "c0"     = 2,                      # prior of baseline hazard
  "tau"    = 0.0375,                 # sd (spike) for coefficient prior
  "cb"     = 20,                     # sd (slab) for coefficient prior
  "pi.ga"  = 0.02,                   # prior variable selection probability for standard Cox models
  "a"      = -4,                     # hyperparameter in MRF prior
  "b"      = 0.1,                    # hyperparameter in MRF prior
  "G"      = simData$G               # hyperparameter in MRF prior
)

# Run a 'Pooled' Bayesian Cox model with graphical learning

set.seed(715074)
BayesSurvive_wrap <- function(use_cpp = FALSE) {
  # suppressWarnings(
    BayesSurvive(
      survObj = dataset, model.type = "Pooled", MRF.G = TRUE, verbose = TRUE,
      hyperpar = hyperparPooled, initial = initial, nIter = 100, burnin = 100,
      cpp = use_cpp
    )
  # )
}
fit_R <- BayesSurvive_wrap(use_cpp = FALSE)
fit_C <- BayesSurvive_wrap(use_cpp = TRUE)

test_that("R and C++ objects are similar", {
  expect_equal(fit_R$call, fit_C$call)
  expect_equal(fit_R$input, fit_C$input)
  for (obj in names(fit_R$output)[2]) {
    expect_equal(fit_R$output[[obj]], fit_C$output[[obj]], tolerance = 1)
  }
})
