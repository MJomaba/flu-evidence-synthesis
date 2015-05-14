#ifndef STATE_HH 
#define STATE_HH

#include<string>
#include<vector>

#include <gsl/gsl_rng.h>

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
    typedef struct {
        int m_plus [260];
        double epsilon[5];
        double psi;
        double transmissibility;
        double susceptibility[7];
        double init_pop;
    } parameter_set ;

/// Keeps the current state of the model/mcmc
struct state_t 
{
    parameter_set parameters;
    double time_infectious, time_latent;

    double positivity_ij[260];
    std::vector<size_t> number_contacts;
};

/// Load an (initial) state from a file
state_t load_state( const std::string &file_path,
        const size_t number_age_groups, const size_t dim_poly_part );

    void save_state(const std::string &file_path, const size_t k, const state_t &current_state, const int * curr_cnt_number, const double *current_contact_regular, const double * result_by_week, const double lv, const double Accept_rate);

};

#endif
