// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Fit
List Fit(arma::mat x0, arma::vec y, arma::vec tV, double lam1, double alpha, double lam2, arma::vec beta, arma::vec Intercept, arma::vec w, arma::vec IndFor0, arma::vec IndTFor0);
RcppExport SEXP _TimePointsPenalized_Fit(SEXP x0SEXP, SEXP ySEXP, SEXP tVSEXP, SEXP lam1SEXP, SEXP alphaSEXP, SEXP lam2SEXP, SEXP betaSEXP, SEXP InterceptSEXP, SEXP wSEXP, SEXP IndFor0SEXP, SEXP IndTFor0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tV(tVSEXP);
    Rcpp::traits::input_parameter< double >::type lam1(lam1SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type lam2(lam2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Intercept(InterceptSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type IndFor0(IndFor0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type IndTFor0(IndTFor0SEXP);
    rcpp_result_gen = Rcpp::wrap(Fit(x0, y, tV, lam1, alpha, lam2, beta, Intercept, w, IndFor0, IndTFor0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TimePointsPenalized_Fit", (DL_FUNC) &_TimePointsPenalized_Fit, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_TimePointsPenalized(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
