#ifndef STATE_HH 
#define STATE_HH 
#include<string>
#include<vector>

#include <gsl/gsl_rng.h>

#include "bson/bson_stream.hh"

#define NAG 7
#define NAG2 49
#define length 364
#define length_weeks 52
#define twopi 6.283185
#define POLY_PART 597
#define h_step 0.25              /*integration step for the ODE system N.B. must be 1/integer and 1/h should be a multiple of 2*/
#define d_app 2
#define dim_par 9
#define dim_par2 81
#define seed 578

/*declaring the random number*/
extern gsl_rng * r;

namespace flu
{

    /// Parameters of the model
    // TODO: Add as vector option?
    struct parameter_set {
        std::vector<double> epsilon = std::vector<double>(5);
        double psi;
        double transmissibility;
        std::vector<double> susceptibility = std::vector<double>(7);
        double init_pop;

        /// Output parameter_set as BSON (JSON)
        friend mongo::BSONEmitter &operator<<(
                mongo::BSONEmitter &bbuild, const parameter_set &pars );

        /// Turn BSON into parameter_set
        friend void operator>>( const mongo::BSONElement &el,
                parameter_set &pars );
    };

    /// Keeps the current state of the model/mcmc
    struct state_t 
    {
        parameter_set parameters;
        double time_infectious, time_latent;

        std::vector<double> positivity_ij = std::vector<double>(260);
        std::vector<size_t> contact_ids;

        /// Output state_t as BSON (JSON)
        friend mongo::BSONEmitter &operator<<(
                mongo::BSONEmitter &bbuild, const state_t &state );

        /// Convert BSON into state_t
        friend void operator>>( const mongo::BSONElement &el,
                state_t &state );
    };

    /// Load an (initial) state from a file
    state_t load_state_json( const std::string &file_path );

    /// Save state to a file
    void save_state_json( const state_t &state, 
            const std::string &file_path );

    state_t load_state( const std::string &file_path, 
        const size_t number_age_groups, const size_t dim_poly_part );

    void save_state(const std::string &file_path, const size_t k, const state_t &state, const std::vector<double> &contact_mat, const double * result_by_week, const double lv, const double Accept_rate);

};

#endif
