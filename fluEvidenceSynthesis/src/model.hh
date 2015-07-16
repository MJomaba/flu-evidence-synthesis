#ifndef MODEL_HH
#define MODEL_HH

#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_NO_DEBUG
#define EIGEN_FAST_MATH

#include<iostream>
#include<cmath>

#include<vector>

#include <boost/date_time.hpp>

#include "state.hh"
#include "vaccine.hh"

#include "rcppwrap.hh"
#include<RcppEigen.h>

namespace flu
{
    boost::posix_time::ptime getTimeFromWeekYear( int week, int year );

    void one_year_SEIR_with_vaccination(double *, const std::vector<double> &pop_vec, double *, const double, const double,  const std::vector<double> &, const Eigen::MatrixXd &contact_regular, double, 
            const vaccine::vaccine_t &vaccine_programme );

    void one_year_SEIR_without_vaccination(double *, const std::vector<double> &pop_vec, double *, double, double, const std::vector<double> &, const Eigen::MatrixXd &contact_regular, double);

    void days_to_weeks(double *, double *);
    void days_to_weeks_no_class(double *, double *);
    void days_to_weeks_5AG(double *, double *);
    double log_likelihood_hyper_poisson(const std::vector<double> &, 
            double, double *,
            Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
            Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
            double *, int);

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
