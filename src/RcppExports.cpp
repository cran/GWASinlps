// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// arma_cor
arma::vec arma_cor(arma::mat X, arma::vec y);
RcppExport SEXP _GWASinlps_arma_cor(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(arma_cor(X, y));
    return rcpp_result_gen;
END_RCPP
}
// looprun
int looprun(CharacterVector varsselected, CharacterVector varsleft, int max_nocollect, double m, int nskip);
RcppExport SEXP _GWASinlps_looprun(SEXP varsselectedSEXP, SEXP varsleftSEXP, SEXP max_nocollectSEXP, SEXP mSEXP, SEXP nskipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type varsselected(varsselectedSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type varsleft(varsleftSEXP);
    Rcpp::traits::input_parameter< int >::type max_nocollect(max_nocollectSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type nskip(nskipSEXP);
    rcpp_result_gen = Rcpp::wrap(looprun(varsselected, varsleft, max_nocollect, m, nskip));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GWASinlps_arma_cor", (DL_FUNC) &_GWASinlps_arma_cor, 2},
    {"_GWASinlps_looprun", (DL_FUNC) &_GWASinlps_looprun, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_GWASinlps(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
