#ifndef STATE_HH 
#define STATE_HH 
#include<string>
#include<vector>

#include<Eigen/Core>

#define POLY_PART 597
//#define seed 578

namespace flu
{
    /// Parameters of the model
    // TODO: Add as vector option?
    struct parameter_set {
        Eigen::VectorXd epsilon;
        double psi;
        double transmissibility;
        Eigen::VectorXd susceptibility;
        double init_pop;
    };

    /// Keeps the current state of the model/mcmc
    struct state_t 
    {
        parameter_set parameters;
        double time_infectious, time_latent;

        std::vector<size_t> contact_ids;
        double likelihood;
    };
}

#endif
