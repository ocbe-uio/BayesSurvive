#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _BayesSurvive_calJpost_helper_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BayesSurvive_matProdVec(SEXP, SEXP);
extern SEXP _BayesSurvive_settingInterval_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _BayesSurvive_sumMatProdVec(SEXP, SEXP);
extern SEXP _BayesSurvive_updateBH_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BayesSurvive_updateBH_list_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BayesSurvive_updateRP_genomic_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BayesSurvive_func_MCMC_graph_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_BayesSurvive_calJpost_helper_cpp",  (DL_FUNC) &_BayesSurvive_calJpost_helper_cpp,   9},
  {"_BayesSurvive_matProdVec",           (DL_FUNC) &_BayesSurvive_matProdVec,            2},
  {"_BayesSurvive_settingInterval_cpp",  (DL_FUNC) &_BayesSurvive_settingInterval_cpp,   4},
  {"_BayesSurvive_sumMatProdVec",        (DL_FUNC) &_BayesSurvive_sumMatProdVec,         2},
  {"_BayesSurvive_updateBH_cpp",         (DL_FUNC) &_BayesSurvive_updateBH_cpp,          7},
  {"_BayesSurvive_updateBH_list_cpp",    (DL_FUNC) &_BayesSurvive_updateBH_list_cpp,     7},
  {"_BayesSurvive_updateRP_genomic_cpp", (DL_FUNC) &_BayesSurvive_updateRP_genomic_cpp, 12},
  {"_BayesSurvive_func_MCMC_graph_cpp",  (DL_FUNC) &_BayesSurvive_func_MCMC_graph_cpp,   6},
  {NULL, NULL, 0}
};

void R_init_BayesSurvive(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
