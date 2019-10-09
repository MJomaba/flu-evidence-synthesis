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
            const Eigen::VectorXd &pop_vec, 
            const Eigen::VectorXd &initial_infected, 
            const double, const double,  
            const Eigen::VectorXd &, 
            const Eigen::MatrixXd &contact_regular, double, 
            const vaccine::vaccine_t &vaccine_programme,
            size_t minimal_resolution = 24,
            const boost::posix_time::ptime &starting_time = 
                getTimeFromWeekYear( 35, 1970 ) );

    /**
     * \brief Run the model for a year
     *
     * @pop_vec Total population in all age groups and risk groups 
     * @initial_infected Total infected population in all age groups and risk groups 
     * @minimal_resolution gives the time resolution (precision) of the returned simulation (in hours)
     */
    cases_t infectionODE(
            const Eigen::VectorXd &pop_vec, 
            const Eigen::VectorXd &initial_infected, 
            const double, const double,  
            const Eigen::VectorXd &, 
            const Eigen::MatrixXd &contact_regular, double, 
            const vaccine::vaccine_t &vaccine_programme,
            size_t minimal_resolution = 24,
            const boost::posix_time::ptime &starting_time = 
                getTimeFromWeekYear( 35, 1970 ) );

    cases_t infectionODE(
            const Eigen::VectorXd &Npop,  
            const Eigen::VectorXd &seed_vec, 
            const double tlatent, const double tinfectious, 
            const Eigen::VectorXd &s_profile, 
            const Eigen::MatrixXd &contact_regular, 
            double transmissibility,
            const vaccine::vaccine_t &vaccine_programme,
            const std::vector<boost::posix_time::ptime> &times );

    void days_to_weeks(double *, double *);
    void days_to_weeks_no_class(double *, double *);

    Eigen::MatrixXd days_to_weeks_AG(const cases_t &simulation,
        const Eigen::MatrixXd &mapping, size_t no_data);

    /// Returns (simplified) log likelihood of one prediction
    long double binomial_log_likelihood( double epsilon, 
            size_t predicted, double population_size, 
            int ili_cases, int ili_monitored,
            int confirmed_positive, int confirmed_samples);

    /// Returns log likelihood of one prediction
    long double log_likelihood(double epsilon, double psi, size_t predicted,
                               double population_size, int ili_cases,
                               int ili_monitored, int confirmed_positive,
                               int confirmed_samples, double abs_err = 1e-5);

    double log_likelihood_hyper_poisson(const Eigen::VectorXd &eps, double psi,
                                        const Eigen::MatrixXd &result_by_week,
                                        const Eigen::MatrixXi &ili,
                                        const Eigen::MatrixXi &mon_pop,
                                        const Eigen::MatrixXi &n_pos,
                                        const Eigen::MatrixXi &n_samples,
                                        const Eigen::VectorXd &pop_AG_RCGP,
                                        double abs_err);

    /**
     * \brief Return the log prior probability of the proposed parameters - current parameters
     *
     * \param susceptibility whether to use the prior based on 2003/04
     */
    double log_prior( const parameter_set &proposed,
            const parameter_set &current,
            bool susceptibility = false );

    /**
     * \brief Return the log prior probability of the proposed parameters - current parameters
     *
     * \param susceptibility whether to use the prior based on 2003/04
     */
    double log_prior( const Eigen::VectorXd &proposed,
            const Eigen::VectorXd &current,
            bool susceptibility = false );
}

#endif
