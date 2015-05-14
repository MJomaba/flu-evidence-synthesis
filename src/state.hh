#ifndef STATE_HH 
#define STATE_HH

#include<string>

#include<model.hh>

namespace flu {

/// Keeps the current state of the model/mcmc
struct state_t 
{
    parameter_set parameters;
    double time_infectious, time_latent;
};

/// Load an (initial) state from a file
state_t load_state( const std::string &file_path,
        const size_t number_age_groups );

};

#endif
