% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.BayesSurvive.R
\name{plot.BayesSurvive}
\alias{plot.BayesSurvive}
\title{Create a plot of estimated coefficients}
\usage{
\method{plot}{BayesSurvive}(x, type = "mean", interval = TRUE, subgroup = 1, ...)
}
\arguments{
\item{x}{an object of class \code{BayesSurvive} or a matrix. If \code{x}
is a matrix, use \code{BayesSurvive:::plot.BayesSurvive(x)}}

\item{type}{type of point estimates of regression coefficients. One of
\code{c("mean", "median")}. Default is \code{mean}}

\item{interval}{logical argument to show 95\% credible intervals. Default
is \code{TRUE}}

\item{subgroup}{index of the subgroup for visualizing posterior coefficients}

\item{...}{additional arguments sent to \code{ggplot2::geom_point()}}
}
\value{
ggplot object
}
\description{
Plot point estimates of regression coefficients and 95\% credible intervals
}
\examples{

library("BayesSurvive")
set.seed(123)

# Load the example dataset
data("simData", package = "BayesSurvive")

dataset <- list(
  "X" = simData[[1]]$X,
  "t" = simData[[1]]$time,
  "di" = simData[[1]]$status
)

# Initial value: null model without covariates
initial <- list("gamma.ini" = rep(0, ncol(dataset$X)))
# Hyperparameters
hyperparPooled <- list(
  "c0"     = 2, # prior of baseline hazard
  "tau"    = 0.0375, # sd for coefficient prior
  "cb"     = 20, # sd for coefficient prior
  "pi.ga"  = 0.02, # prior variable selection probability for standard Cox models
  "a"      = -4, # hyperparameter in MRF prior
  "b"      = 0.1, # hyperparameter in MRF prior
  "G"      = simData$G # hyperparameter in MRF prior
)

\donttest{
# run Bayesian Cox with graph-structured priors
fit <- BayesSurvive(
  survObj = dataset, hyperpar = hyperparPooled,
  initial = initial, nIter = 50
)

# show posterior mean of coefficients and 95\% credible intervals
library("GGally")
plot(fit) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, size = 7))
}

}
