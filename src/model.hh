#ifndef MODEL_HH
#define MODEL_HH


#include<iostream>
#include<cmath>

#include<vector>

#include "state.hh"
#include "io.hh"

namespace flu
{
    void one_year_SEIR_with_vaccination(double *, const std::vector<double> &pop_vec, double *, const double, const double, double const*, double *, double, double *,double *);
    void one_year_SEIR_without_vaccination(double *, const std::vector<double> &pop_vec, double *, double, double, double const*, double *, double);
    void days_to_weeks(double *, double *);
    void days_to_weeks_no_class(double *, double *);
    void days_to_weeks_5AG(double *, double *);
    double log_likelihood_b(double *, double *, double *, double *, int *, int *, double *);
    double log_likelihood_hyper(double *, double *, int *, int *, int *, int *, double *);
    double log_likelihood_hyper_poisson(double *, double, double *, int *, int *, int *, int *, double *, int);

    void save_scenarii( FILE *Scen1FS, FILE *Scen2FS, const std::vector<double> &pop_vec,  double *prop_init_inf, const state_t &state, double *contact_mat, int n_scenarii, double **vaccine_cal, double **vaccine_efficacy_year, std::string path, int * First_write );

    void proposal_haario(parameter_set *, parameter_set *, double *, double *, int , double);
    void proposal_haario_adapt_scale(parameter_set *, parameter_set *, double *, double *, int n, double, double);
    void cholevsky(double *, double *, int);
    void update_sum_corr(double *, parameter_set *);
};

#endif
