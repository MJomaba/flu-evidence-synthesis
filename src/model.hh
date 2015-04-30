#ifndef MODEL_HH
#define MODEL_HH

#include <gsl/gsl_rng.h>

#include<iostream>
#include<cmath>

#include "boost/filesystem.hpp"

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


typedef struct {
    int m_plus [260];
    double epsilon[5];
    double psi;
    double transmissibility;
    double susceptibility[7];
    double init_pop;
} parameter_set ;

void one_year_SEIR_with_vaccination(double *, double *, double *, double, double, double *, double *, double, double *,double *);
void one_year_SEIR_without_vaccination(double *, double *, double *, double, double, double *, double *, double);
void days_to_weeks(double *, double *);
void days_to_weeks_no_class(double *, double *);
void days_to_weeks_5AG(double *, double *);
double log_likelihood_b(double *, double *, double *, double *, int *, int *, double *);
double log_likelihood_hyper(double *, double *, int *, int *, int *, int *, double *);
double log_likelihood_hyper_poisson(double *, double, double *, int *, int *, int *, int *, double *, int);
double normal(double, double);
void save_state(const char *, int, double, double, double, double, double *, double *, double *, double, int *, double *, double *, double, double);
void save_scenarii( FILE *, FILE *, int, double *,  double *,  double, double, double, double *, double *, int, double **, double **, char *, int *);
void proposal_haario(parameter_set *, parameter_set *, double *, double *, int , double);
void proposal_haario_adapt_scale(parameter_set *, parameter_set *, double *, double *, int n, double, double);
void cholevsky(double *, double *, int);
void update_sum_corr(double *, parameter_set *);

FILE * read_file( const std::string path, const std::string filename );
FILE * write_file( const std::string filename );
FILE * append_file( const std::string filename );

#endif
