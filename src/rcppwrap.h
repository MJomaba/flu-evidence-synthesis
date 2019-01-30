#ifndef RCPP_WRAP_HH
#define RCPP_WRAP_HH

#include <RcppCommon.h>
#include "vaccine.h"
#include "state.h"
#include "contacts.h"
#include "inference.h"

namespace Rcpp {
    using namespace flu;

    template <> parameter_set as( SEXP );
    template <> SEXP wrap( const parameter_set &parameters );
    
    template <> state_t as( SEXP );
    template <> SEXP wrap( const state_t &sample );
    /// Keeps the current state of the model/mcmc
    /*struct state_t 
    {
        parameter_set parameters;
        double time_infectious, time_latent;

        std::vector<double> positivity_ij = std::vector<double>(260);
        std::vector<size_t> contact_ids;
    };*/

    using namespace flu::vaccine;
    template <> vaccine_t as( SEXP );

    using namespace flu::contacts;
    template <> contacts_t as( SEXP );

    template <> SEXP wrap( const mcmc_result_inference_t &mcmcResult );
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
//// [[Rcpp::depends(RcppBDT)]]
#include <RcppEigen.h>
//#include <RcppBDT.h>
#include<Rcpp.h>
#endif
