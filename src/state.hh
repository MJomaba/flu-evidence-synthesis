#ifndef STATE_HH 
#define STATE_HH

#include<string>
#include<vector>

#include<model.hh>

namespace flu {

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
