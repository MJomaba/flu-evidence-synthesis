#ifndef RCPP_WRAP_HH
#define RCPP_WRAP_HH

#include <RcppCommon.h>
#include "vaccine.hh"
#include "state.hh"

using namespace flu;
using namespace flu::vaccine;
namespace Rcpp {
    template <> vaccine_t as( SEXP );

    template <> parameter_set as( SEXP );
    /* 
    /// Parameters of the model
    struct parameter_set {
        std::vector<double> epsilon = std::vector<double>(5);
        double psi;
        double transmissibility;
        std::vector<double> susceptibility = std::vector<double>(7);
        double init_pop;
    };*/

    template <> state_t as( SEXP );
    /// Keeps the current state of the model/mcmc
    /*struct state_t 
    {
        parameter_set parameters;
        double time_infectious, time_latent;

        std::vector<double> positivity_ij = std::vector<double>(260);
        std::vector<size_t> contact_ids;
    };*/
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include<Rcpp.h>
#endif
