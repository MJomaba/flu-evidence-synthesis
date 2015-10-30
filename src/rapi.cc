#include <boost/date_time.hpp>

#include "rcppwrap.h"

#include "mcmc.h"

#include "proposal.h"
#include "contacts.h"
#include "model.h"
#include "ode.h"
#include "inference.h"

namespace bt = boost::posix_time;

/**
 * This file contains some thin wrappers around C++ code. Mostly used
 * for testing the C++ code from R.
 *
 * The wrappers are needed, because in the C++ code we often want to pass
 * by const reference, which is impossible in a function called from R
 */

//' Update means when a new posterior sample is calculated
//'
//' @param means the current means of the parameters
//' @param v the new parameter values
//' @param n The number of posterior (mcmc) samples taken till now
//' @return The updated means given the new parameter sample
//'
// [[Rcpp::export(name=".updateMeans")]]
Eigen::VectorXd updateMeans( Eigen::VectorXd means,
        Eigen::VectorXd v, size_t n )
{
    return flu::proposal::updateMeans( means, v, n );
}

//' Update covariance matrix of posterior parameters
//'
//' Used to enable faster mixing of the mcmc chain
//' @param cov The current covariance matrix
//' @param v the new parameter values
//' @param means the current means of the parameters
//' @param n The number of posterior (mcmc) samples taken till now
//' @return The updated covariance matrix given the new parameter sample
//'
// [[Rcpp::export(name=".updateCovariance")]]
Eigen::MatrixXd updateCovariance( Eigen::MatrixXd cov, 
        Eigen::VectorXd v, Eigen::VectorXd means, size_t n )
{
    return flu::proposal::updateCovariance( cov, v, means, n );
}

//' Convert given week in given year into an exact date corresponding to the Monday of that week
//'
//' @param week The number of the week we need the date of
//' @param year The year
//' @return The date of the Monday in that week 
//'
// [[Rcpp::export]]
Rcpp::Datetime getTimeFromWeekYear( int week, int year )
{
    //return bt::to_iso_string(flu::getTimeFromWeekYear( week, year ));
    return Rcpp::Datetime(
            bt::to_iso_extended_string(flu::getTimeFromWeekYear( week, year )),
            "%Y-%m-%dT%H:%M:%OS");
    //return Rcpp::wrap(flu::getTimeFromWeekYear( week, year ));
    //return Rcpp::wrap(flu::getTimeFromWeekYear( week, year ));
}

//' Run the SEIR model for the given parameters
//'
//' @param age_sizes A vector with the population size by each age {1,2,..}
//' @param vaccine_calendar A vaccine calendar valid for that year
//' @param polymod_data Contact data for different age groups
//' @param susceptibility Vector with susceptibilities of each age group
//' @param transmissibility The transmissibility of the strain
//' @param init_pop The (log of) initial infected population
//' @param infection_delays Vector with the time of latent infection and time infectious
//' @param interval Interval (in days) between data points
//' @return A data frame with number of new cases after each interval during the year
//'
// [[Rcpp::export(name="infection.model")]]
Rcpp::DataFrame runSEIRModel(
        std::vector<size_t> age_sizes, 
        flu::vaccine::vaccine_t vaccine_calendar,
        flu::contacts::contacts_t polymod_data,
        Eigen::VectorXd susceptibility, 
        double transmissibility, 
        double init_pop,
        Eigen::VectorXd infection_delays, 
        size_t interval = 1 )
{

    flu::data::age_data_t age_data;
    age_data.age_sizes = age_sizes;
    age_data.age_group_sizes = flu::data::group_age_data( age_sizes );
    auto pop_vec = flu::data::separate_into_risk_groups( 
            age_data.age_group_sizes );

    /*translate into an initial infected population*/
    auto curr_init_inf = Eigen::VectorXd::Constant( 
            susceptibility.size(), pow(10,init_pop) );

    auto current_contact_regular = 
        flu::contacts::to_symmetric_matrix( polymod_data, age_data );

    auto result = flu::one_year_SEIR_with_vaccination(pop_vec, curr_init_inf, infection_delays[0], infection_delays[1], susceptibility, current_contact_regular, transmissibility, vaccine_calendar, interval*24 );

    Rcpp::List resultList( result.cases.cols() + 1 );
    Rcpp::CharacterVector columnNames;

    //Rcpp::DataFrame densities = Rcpp::wrap<Rcpp::DataFrame>( result.cases );
    // Convert times
    auto times = Rcpp::DatetimeVector(result.times.size());
    for ( size_t i = 0; i < result.times.size(); ++i )
    {
        times[i] = 
                Rcpp::Datetime(
            bt::to_iso_extended_string( result.times[i] ),
            "%Y-%m-%dT%H:%M:%OS");
    }

    columnNames.push_back( "Time" );
    resultList[0] = times;

    for (int i=0; i<result.cases.cols(); ++i)
    {
        resultList[i+1] = Eigen::VectorXd(result.cases.col(i));
        columnNames.push_back( 
                "V" + boost::lexical_cast<std::string>( i+1 ) );
    }

    resultList.attr("names") = columnNames;

    return Rcpp::DataFrame(resultList);
    //return densities;
    //return resultMatrix;
}


//' Returns log likelihood of the predicted number of cases given the data for that week
//'
//' The model results in a prediction for the number of new cases in a certain age group and for a certain week. This function calculates the likelihood of that given the data on reported Influenza Like Illnesses and confirmed samples.
//'
//' @param epsilon Parameter for the probability distribution
//' @param psi Parameter for the probability distribution
//' @param predicted Number of cases predicted by your model
//' @param population_size The total population size in the relevant age group
//' @param ili_cases The number of Influenza Like Illness cases
//' @param ili_monitored The size of the population monitored for ILI
//' @param confirmed_positive The number of samples positive for the Influenza strain
//' @param confirmed_samples Number of samples tested for the Influenza strain
//'
// [[Rcpp::export(name="llikelihood.cases")]]
double log_likelihood( double epsilon, double psi, 
        size_t predicted, double population_size, 
        int ili_cases, int ili_monitored,
        int confirmed_positive, int confirmed_samples )
{
    return flu::log_likelihood( epsilon, psi,
            predicted, population_size,
            ili_cases, ili_monitored,
            confirmed_positive, confirmed_samples, 2 );
}


//' Run an ODE model with the runge-kutta solver for testing purposes
//'
//' @param step_size The size of the step between returned time points
//' @param h_step The starting integration delta size
//'
// [[Rcpp::export(name=".runRKF")]]
Eigen::MatrixXd runPredatorPrey(double step_size = 0.1, double h_step=0.01)
{
    Eigen::MatrixXd result(201,3);
    double ct = 0;
    Eigen::VectorXd y(2); y[0] = 10.0; y[1] = 5.0;
    size_t row_count = 0;
    while( ct <= 20.01 )
    {
        result( row_count, 0 ) = ct;
        result( row_count, 1 ) = y[0];
        result( row_count, 2 ) = y[1];
        double t = 0;
        while (t < step_size)
        {
            auto ode_func = []( const Eigen::VectorXd &y, const double ct )
            {
                Eigen::VectorXd dy(2);
                dy[0] = 1.5*y[0]-1.0*y[0]*y[1];
                dy[1] = 1.0*y[0]*y[1]-3.0*y[1];
                return dy;
            };
            y = ode::rkf45_astep( std::move(y), ode_func,
                        h_step, t, step_size, 1.0e-12 );
        }
        ct += step_size;
        ++row_count;
    }
    return result;
}

//' Run an ODE model with the simple step wise solver for testing purposes
//'
//' @param step_size The size of the step between returned time points
//' @param h_step The starting integration delta size
//'
// [[Rcpp::export(name=".runStep")]]
Eigen::MatrixXd runPredatorPreySimple(double step_size = 0.1, double h_step=1e-5)
{
    Eigen::MatrixXd result(201,3);
    double ct = 0;
    Eigen::VectorXd y(2); y[0] = 10.0; y[1] = 5.0;
    size_t row_count = 0;
    while( ct <= 20.01 )
    {
        result( row_count, 0 ) = ct;
        result( row_count, 1 ) = y[0];
        result( row_count, 2 ) = y[1];
        double t = 0;
        while (t < step_size)
        {
            auto ode_func = []( const Eigen::VectorXd &y, const double ct )
            {
                Eigen::VectorXd dy(2);
                dy[0] = 1.5*y[0]-1.0*y[0]*y[1];
                dy[1] = 1.0*y[0]*y[1]-3.0*y[1];
                return dy;
            };
            y = ode::step( std::move(y), ode_func,
                        h_step, t, step_size );
        }
        ct += step_size;
        ++row_count;
    }
    return result;
}

//' Adaptive MCMC algorithm implemented in C++
//'
//' MCMC which adapts its proposal distribution for faster convergence following:
//' Sherlock, C., Fearnhead, P. and Roberts, G.O. The Random Walk Metrolopois: Linking Theory and Practice Through a Case Study. Statistical Science 25, no.2 (2010): 172-190.
//'
//' @param lprior A function returning the log prior probability of the parameters 
//' @param llikelihood A function returning the log likelihood of the parameters given the data
//' @param nburn Number of iterations of burn in
//' @param initial Vector with starting parameter values
//' @param nbatch Number of batches to run (number of samples to return)
//' @param blen Length of each batch
//' 
//' @return Returns a list with the accepted samples and the corresponding llikelihood values
//'
//' @seealso \code{\link{adaptive.mcmc}} For a more flexible R frontend to this function.
//'
// [[Rcpp::export(name="adaptive.mcmc.cpp")]]
Rcpp::List adaptiveMCMCR( 
        Rcpp::Function lprior, Rcpp::Function llikelihood,
        size_t nburn,
        Eigen::VectorXd initial, 
        size_t nbatch, size_t blen = 1 )
{
    auto cppLprior = [&lprior]( const Eigen::VectorXd &pars ) {
        double lPrior = Rcpp::as<double>(lprior( pars ));
        return lPrior;
    };

    auto cppLlikelihood = [&llikelihood]( const Eigen::VectorXd &pars ) {
        double ll = Rcpp::as<double>(llikelihood( pars ));
        return ll;
    };

    auto mcmcResult = flu::adaptiveMCMC( cppLprior, cppLlikelihood, nburn, initial, nbatch, blen );
    Rcpp::List rState;
    rState["batch"] = Rcpp::wrap( mcmcResult.batch );
    rState["llikelihoods"] = Rcpp::wrap( mcmcResult.llikelihoods );
    return rState;
}
