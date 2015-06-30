#ifndef RCPP_WRAP_HH
#define RCPP_WRAP_HH

#include <RcppCommon.h>
#include "vaccine.hh"
#include "state.hh"
#include "contacts.hh"

namespace Rcpp {
    using namespace flu;

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

    using namespace flu::vaccine;
    template <> vaccine_t as( SEXP );

    using namespace flu::contacts;
    template <> contacts_t as( SEXP );
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
#include<Rcpp.h>
#endif
