[![R-CMD-check](https://github.com/ocbe-uio/BayesSurvive/workflows/R-CMD-check/badge.svg)](https://github.com/ocbe-uio/BayesSurv/actions)
[![License](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)


# BayesSurvive


This is a R/Rcpp package **BayesSurvive** for Bayesian survival models with graph-structured selection priors for sparse identification of high-dimensional features predictive of survival ([Madjar et al., 2021](https://doi.org/10.1186/s12859‐021‐04483‐z)) and its extensions with the use of a fixed graph via a Markov Random Field (MRF) prior for capturing known structure of high-dimensional features, e.g. disease-specific pathways from the Kyoto Encyclopedia of Genes and Genomes (KEGG) database.

## Installation

Install the latest development version from GitHub

```r
#install.packages("remotes")
remotes::install_github("ocbe-uio/BayesSurvive")
```

## Examples

### Simulate data

```r
library("BayesSurvive")
# Load the example dataset
data("simData", package = "BayesSurvive")
dataset = list("X" = simData[[1]]$X, 
               "t" = simData[[1]]$time,
               "di" = simData[[1]]$status)
```

### Run a Bayesian Cox model

```r
## Initial value: null model without covariates
initial = list("gamma.ini" = rep(0, ncol(dataset$X)))
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

## run Bayesian Cox with graph-structured priors
fit <- BayesSurvive(survObj = dataset, model.type = "Pooled", MRF.G = TRUE, 
                    hyperpar = hyperparPooled, initial = initial, nIter = 100)

## show posterior mean of coefficients and 95% credible intervals
library("GGally")
plot(fit) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 7))

#plot(fit$output$beta.p[,1], type="l")
#fit$output$beta.margin
#fit$output$gamma.margin
#simData[[1]]$trueB
```

<img src="man/figures/README_plot_beta.png" width="100%" />


### Plot time-dependent Brier scores

The function `BayesSurvive::plotBrier()` can show the time-dependent Brier scores based on posterior mean of coefficients or Bayesian model averaging.

```r
plotBrier(fit, , survObj.new = dataset)
```

<img src="man/figures/README_plot_brier.png" width="50%" />

The integrated Brier score (IBS) can be obtained by the function `BayesSurvive::predict()`.

```r
predict(fit, survObj.new = dataset)
```
```{ .text .no-copy }
    Null.model Bayesian.Cox
IBS 0.09147208   0.03797003
```

### Predict survival probabilities and cumulative hazards

The function `BayesSurvive::predict()` can estimate the survival probabilities and cumulative hazards.

```r
predict(fit, survObj.new = dataset, type = c("cumhazard", "survival"))
```
```{ .text .no-copy }
       observation times cumhazard survival
    1:           1   3.3  1.82e-04 1.00e+00
    2:           2   3.3  2.42e-01 7.85e-01
    3:           3   3.3  2.93e-06 1.00e+00
    4:           4   3.3  5.40e-03 9.95e-01
    5:           5   3.3  8.52e-04 9.99e-01
   ---                                     
 9996:          96   9.5  1.38e+01 9.68e-07
 9997:          97   9.5  9.71e+01 6.94e-43
 9998:          98   9.5  1.35e+00 2.58e-01
 9999:          99   9.5  2.24e+00 1.06e-01
10000:         100   9.5  8.62e+00 1.81e-04
```

### Run a 'Pooled' Bayesian Cox model with graphical learning

```r
hyperparPooled <- append(hyperparPooled, list("lambda" = 3, "nu0" = 0.05, "nu1" = 5))
fit2 <- BayesSurvive(survObj = list(dataset), model.type = "Pooled", MRF.G = FALSE,
                     hyperpar = hyperparPooled, initial = initial, nIter = 10)
```

### Run a Bayesian Cox model with subgroups using fixed graph 

```r
# specify a fixed joint graph between two subgroups
hyperparPooled$G <- Matrix::bdiag(simData$G, simData$G)
dataset2 <- simData[1:2]
dataset2 <- lapply(dataset2, setNames, c("X", "t", "di", "X.unsc", "trueB"))
fit3 <- BayesSurvive(survObj = dataset2, 
                     hyperpar = hyperparPooled, initial = initial, 
                     model.type="CoxBVSSL", MRF.G = TRUE, 
                     nIter = 10, burnin = 5)
```

### Run a Bayesian Cox model with subgroups using graphical learning

```r
fit4 <- BayesSurvive(survObj = dataset2, 
                     hyperpar = hyperparPooled, initial = initial, 
                     model.type="CoxBVSSL", MRF.G = FALSE, 
                     nIter = 3, burnin = 0)
```

## References

> Katrin Madjar, Manuela Zucknick, Katja Ickstadt, Jörg Rahnenführer (2021).
> Combining heterogeneous subgroups with graph‐structured variable selection priors for Cox regression.
> _BMC Bioinformatics_, 22(1):586. DOI:[10.1186/s12859-021-04483-z](https://doi.org/10.1186/s12859-021-04483-z).
