#ifndef MODEL_HH
#define MODEL_HH


#include<iostream>
#include<cmath>

#include<vector>

#include <boost/numeric/ublas/matrix.hpp>

#include "state.hh"
#include "vaccine.hh"

namespace flu
{
    namespace bu = boost::numeric::ublas;

    void one_year_SEIR_with_vaccination(double *, const std::vector<double> &pop_vec, double *, const double, const double,  const std::vector<double> &, const bu::matrix<double> &contact_regular, double, 
            const vaccine::vaccine_t &vaccine_programme );

    void one_year_SEIR_without_vaccination(double *, const std::vector<double> &pop_vec, double *, double, double, const std::vector<double> &, const bu::matrix<double> &contact_regular, double);

    void days_to_weeks(double *, double *);
    void days_to_weeks_no_class(double *, double *);
    void days_to_weeks_5AG(double *, double *);
    double log_likelihood_hyper_poisson(const std::vector<double> &, double, double *, int *, int *, int *, int *, double *, int);

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
