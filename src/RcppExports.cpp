// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// kendalltaudist
Rcpp::NumericMatrix kendalltaudist(const Rcpp::NumericMatrix& V);
RcppExport SEXP _GSReg_kendalltaudist(SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(kendalltaudist(V));
    return rcpp_result_gen;
END_RCPP
}
// kendalltaudistFromTemp
Rcpp::NumericVector kendalltaudistFromTemp(const Rcpp::NumericMatrix& V, const Rcpp::NumericMatrix& T);
RcppExport SEXP _GSReg_kendalltaudistFromTemp(SEXP VSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(kendalltaudistFromTemp(V, T));
    return rcpp_result_gen;
END_RCPP
}
// kendalltaudistRestricted
Rcpp::NumericMatrix kendalltaudistRestricted(const Rcpp::NumericMatrix& V, const Rcpp::NumericMatrix& R);
RcppExport SEXP _GSReg_kendalltaudistRestricted(SEXP VSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(kendalltaudistRestricted(V, R));
    return rcpp_result_gen;
END_RCPP
}
// Nij
Rcpp::NumericMatrix Nij(const Rcpp::NumericMatrix& V);
RcppExport SEXP _GSReg_Nij(SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(Nij(V));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GSReg_kendalltaudist", (DL_FUNC) &_GSReg_kendalltaudist, 1},
    {"_GSReg_kendalltaudistFromTemp", (DL_FUNC) &_GSReg_kendalltaudistFromTemp, 2},
    {"_GSReg_kendalltaudistRestricted", (DL_FUNC) &_GSReg_kendalltaudistRestricted, 2},
    {"_GSReg_Nij", (DL_FUNC) &_GSReg_Nij, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_GSReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
