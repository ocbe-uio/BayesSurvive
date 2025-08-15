# BayesSurvive (development version)

* Fixed code smells (issue #24)
* Use `calJpost_cpp()` in all `cpp == TRUE` cases (issue #11)
* Fix `lambda` and `pi_G` assignment for `method == "CovBVSSL"` (issue #11)
* Translated `calJpost()` for `method == "Pooled" && MRF_G` (issue #28)

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
