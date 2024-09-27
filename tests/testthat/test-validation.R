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
hyperparPooled <- list(
  "c0"     = 2, # prior of baseline hazard
  "tau"    = 0.0375, # sd (spike) for coefficient prior
  "cb"     = 20, # sd (slab) for coefficient prior
  "pi.ga"  = 0.02, # prior variable selection probability for standard Cox models
  "a"      = -4, # hyperparameter in MRF prior
  "b"      = 0.1, # hyperparameter in MRF prior
  "G"      = simData$G # hyperparameter in MRF prior
)

test_that("burnin must be less than nIter", {
  expect_error(BayesSurvive(dataset, Iter = 30, burnin = 30))
})
