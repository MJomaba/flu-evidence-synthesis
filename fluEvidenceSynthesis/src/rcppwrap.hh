#ifndef RCPP_WRAP_HH
#define RCPP_WRAP_HH

#include <RcppCommon.h>
#include "vaccine.hh"

using namespace flu::vaccine;
namespace Rcpp {
    template <> vaccine_t as( SEXP );
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include<Rcpp.h>
#endif
