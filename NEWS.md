# BayesSurvive (development version)

# BayesSurvive 0.1.0

* Add test code
* Translate R functions `func_MCMC_graph()`, `updateGamma()` and `UpdateRPlee11()` to C++ functions
* Add `misc.cpp` including functions to convert `Rcpp::List` to armadillo vector, matrix or cube
* Update function `VS()` to allow a matrix/array input and allow multiple subgroups
* Rename the output of MPM coefficients in function `coef.BayesSurvive()`
* Add the index of prediction accuracy (IPA) in function `predict.BayesSurvive()`

# BayesSurvive 0.0.2

* Add vignette
* Add function `coef.BayesSurvive()` to extract estimated coefficients and uncertainty
* Add function `VS()` to determine variable selection
* Fix documentation issues

# BayesSurvive 0.0.1

* First version
