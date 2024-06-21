# Load the example dataset
dataset <- list(
  "X" = simData[[1]]$X,
  "t" = simData[[1]]$time,
  "di" = simData[[1]]$status
)

### Run a Bayesian Cox model

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

## run Bayesian Cox with graph-structured priors
set.seed(123)
fit <- BayesSurvive(
  survObj = dataset, model.type = "Pooled", MRF.G = TRUE,
  hyperpar = hyperparPooled, initial = initial,
  nIter = 50, burnin = 0,
  verbose = FALSE
)
pred_1 <- predict(fit, survObj.new = dataset, times = 8.5, verbose = FALSE)
pred_2 <- predict(fit, survObj.new = dataset, type = c("cumhazard", "survival"))

test_that("fit has properly class and length", {
  expect_s3_class(fit, "BayesSurvive")
  expect_length(fit, 3L)
  expect_length(fit$input, 9L)
  expect_length(fit$input$hyperpar, 10L)
  expect_length(fit$output, 17L)
  expect_length(fit$output$survObj, 4L)
})

test_that("fit has expected values", {
  tol <- 1e-3
  with(fit$output, {
    expect_equal(eta0, c("(Intercept)" = 1.74e-5), tolerance = tol)
    expect_equal(head(s, 4), c(3.2969, 3.3217, 4.0938, 4.4107), tolerance = tol)
    expect_equal(head(survObj$t, 4), c(8.53, 4.09, 8.82, 6.09), tolerance = tol)
  })
  expect_equal(which(VS(fit, method = "FDR", threshold = 0.6)), c(6, 12))
})

test_that("predictions have expected values", {
  tol <- 1e-1
  expect_equal(
    head(pred_1$times),
    c(0.00000000, 0.08585859, 0.17171717, 0.25757576, 0.34343434, 0.42929293),
    tolerance = tol
  )
  expect_equal(
    head(pred_2$cumhazard[, 1], 5L),
    c(1.782968e-04, 1.544469e-01, 9.193403e-07, 6.192432e-03, 4.701002e-04),
    tolerance = tol
  )
  expect_equal(
    head(pred_2$survival[, 1]),
    c(0.9998217, 0.8568890, 0.9999991, 0.9938267, 0.9995300, 0.9935771),
    tolerance = tol
  )
  expect_equal(
    head(pred_2$times),
    c(3.296921, 3.321787, 3.899565, 4.093849, 4.410782, 4.932070),
    tolerance = tol
  )
  expect_false(pred_2$se)
  expect_false(pred_2$band)
  expect_false(pred_2$diag)
  expect_false(pred_2$baseline)
})

### Run a 'Pooled' Bayesian Cox model with graphical learning

hyperparPooled <- append(hyperparPooled, list("lambda" = 3, "nu0" = 0.05, "nu1" = 5))
set.seed(3346141)
fit2 <- BayesSurvive(
  survObj = list(dataset), model.type = "Pooled", MRF.G = FALSE,
  hyperpar = hyperparPooled, initial = initial, nIter = 3, verbose = FALSE
)
test_that("fit2 has expected values", {
  tol <- 1e-3
  with(fit2$output, {
    expect_equal(eta0[[1]], c("(Intercept)" = 1.74e-5), tolerance = tol)
    expect_equal(head(s[[1]], 3), c(3.2969, 3.3217, 4.0938), tolerance = tol)
    expect_equal(head(survObj[[1]]$t, 3), c(8.53, 4.09, 8.82), tolerance = tol)
  })
  expect_equal(which(VS(fit2, method = "FDR", threshold = 0.8)), c(81, 182))
})
