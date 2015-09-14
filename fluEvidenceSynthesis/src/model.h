#ifndef MODEL_HH
#define MODEL_HH


#include<iostream>
#include<cmath>

#include<vector>

#include <boost/date_time.hpp>

#include "state.h"
#include "vaccine.h"

#include "rcppwrap.h"
#include<RcppEigen.h>

namespace flu
{
    boost::posix_time::ptime getTimeFromWeekYear( int week, int year );

    struct cases_t 
    {
        //! Number of new cases
        Eigen::MatrixXd cases;

        //! Times corresponding to number of cases
        std::vector<boost::posix_time::ptime> times;
    };

    /**
     * \brief Run the model for a year
     *
     * @minimal_resolution gives the time resolution (precision) of the returned simulation (in hours)
     */
    cases_t one_year_SEIR_with_vaccination(
            const std::vector<double> &pop_vec, 
            double *, const double, const double,  
            const std::vector<double> &, 
            const Eigen::MatrixXd &contact_regular, double, 
            const vaccine::vaccine_t &vaccine_programme,
            size_t minimal_resolution = 24,
            const boost::posix_time::ptime &starting_time = 
                getTimeFromWeekYear( 35, 1970 ) );

    void days_to_weeks(double *, double *);
    void days_to_weeks_no_class(double *, double *);

    Eigen::MatrixXd days_to_weeks_5AG(const cases_t &simulation);

    /// Returns log likelihood of one prediction
    long double log_likelihood( double epsilon, double psi, 
            size_t predicted, double population_size, 
            int ili_cases, int ili_monitored,
            int confirmed_positive, int confirmed_samples, 
            int depth = 2 );

    double log_likelihood_hyper_poisson(const std::vector<double> &eps, 
            double psi, const Eigen::MatrixXd &result_by_week,
            const Eigen::MatrixXi &ili, const Eigen::MatrixXi &mon_pop, 
            const Eigen::MatrixXi &n_pos, const Eigen::MatrixXi &n_samples, 
            double * pop_5AG_RCGP, int depth);

    /**
     * \brief Return the log prior probability of the proposed parameters - current parameters
     *
     * \param susceptibility whether to use the prior based on 2003/04
     */
    double log_prior( const parameter_set &proposed,
            const parameter_set &current,
            bool susceptibility = false );

}

#endif
