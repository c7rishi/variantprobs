// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// good_turing_multinom
Rcpp::List good_turing_multinom(NumericVector Obs, NumericVector Freq, double confid_factor);
RcppExport SEXP _variantprobs_good_turing_multinom(SEXP ObsSEXP, SEXP FreqSEXP, SEXP confid_factorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Obs(ObsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Freq(FreqSEXP);
    Rcpp::traits::input_parameter< double >::type confid_factor(confid_factorSEXP);
    rcpp_result_gen = Rcpp::wrap(good_turing_multinom(Obs, Freq, confid_factor));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_variantprobs_good_turing_multinom", (DL_FUNC) &_variantprobs_good_turing_multinom, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_variantprobs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
