#ifndef MODEL_HH
#define MODEL_HH


#include<iostream>
#include<cmath>

#include<vector>

#include "state.hh"
#include "io.hh"
#include "vaccine.hh"

namespace flu
{
    void one_year_SEIR_with_vaccination(double *, const std::vector<double> &pop_vec, double *, const double, const double, double const*, const std::vector<double> &contact_regular, double, 
            const vaccine::vaccine_t &vaccine_programme );

    void one_year_SEIR_without_vaccination(double *, const std::vector<double> &pop_vec, double *, double, double, double const*, const std::vector<double> &contact_regular, double);

    void days_to_weeks(double *, double *);
    void days_to_weeks_no_class(double *, double *);
    void days_to_weeks_5AG(double *, double *);
    double log_likelihood_b(double *, double *, double *, double *, int *, int *, double *);
    double log_likelihood_hyper(double *, double *, int *, int *, int *, int *, double *);
    double log_likelihood_hyper_poisson(double *, double, double *, int *, int *, int *, int *, double *, int);

    void save_scenarii( FILE *Scen1FS, FILE *Scen2FS, const std::vector<double> &pop_vec,  double *prop_init_inf, const state_t &state, const std::vector<double> &contact_mat, const std::vector<vaccine::vaccine_t> &vaccine_scenarios, std::string path );

    /**
     * \brief Return the log prior probability of the proposed parameters - current parameters
     *
     * \param susceptibility whether to use the prior based on 2003/04
     */
    double log_prior( const parameter_set &proposed,
            const parameter_set &current,
            bool susceptibility = false );

};

#endif
