#include <boost/date_time.hpp>
#include <regex>

#include "rcppwrap.h"

#include "mcmc.h"

#include "proposal.h"
#include "contacts.h"
#include "model.h"
#include "ode.h"
#include "inference.h"
#include "data.h"

namespace bt = boost::posix_time;

Rcpp::Datetime ptime_to_datetime( const bt::ptime &ti )
{
    return Rcpp::Datetime(
            bt::to_iso_extended_string( ti ),
            "%Y-%m-%dT%H:%M:%OS");
}

Rcpp::Date ptime_to_date( const bt::ptime &ti )
{
    return Rcpp::Date(
            bt::to_iso_extended_string( ti ),
            "%Y-%m-%dT%H:%M:%OS");
}


bt::ptime datetime_to_ptime( const Rcpp::Datetime &ti )
{
    return bt::ptime( boost::gregorian::date(ti.getYear(), ti.getMonth(), 
                ti.getDay()),
            bt::hours( ti.getHours() )+
            bt::minutes( ti.getMinutes() ) +
            bt::seconds( ti.getSeconds() )
            );
}

bt::ptime date_to_ptime( const Rcpp::Date &ti )
{
    return bt::ptime( boost::gregorian::date(ti.getYear(), ti.getMonth(), 
                ti.getDay()), bt::hours(12) );
}


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
    return ptime_to_datetime(flu::getTimeFromWeekYear( week, year ));
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
// [[Rcpp::export(name=".infection.model")]]
Rcpp::DataFrame runSEIRModel(
        std::vector<size_t> age_sizes, 
        flu::vaccine::vaccine_t vaccine_calendar,
        Eigen::MatrixXi polymod_data,
        Eigen::VectorXd susceptibility, 
        double transmissibility, 
        double init_pop,
        Eigen::VectorXd infection_delays, 
        size_t interval = 1 )
{

    flu::data::age_data_t age_data;
    age_data.age_sizes = age_sizes;

    std::vector<size_t> age_group_limits = {1,5,15,25,45,65};
    age_data.age_group_sizes = flu::data::group_age_data( age_sizes, 
            age_group_limits );

    Eigen::MatrixXd risk_proportions = Eigen::MatrixXd( 
            2, age_data.age_group_sizes.size() );
    risk_proportions << 
        0.021, 0.055, 0.098, 0.087, 0.092, 0.183, 0.45, 
        0, 0, 0, 0, 0, 0, 0;

    auto pop_vec = flu::data::separate_into_risk_groups( 
            age_data.age_group_sizes, risk_proportions );

    /*translate into an initial infected population*/
    auto curr_init_inf = Eigen::VectorXd::Constant( 
            susceptibility.size(), pow(10,init_pop) );

    auto current_contact_regular = 
        flu::contacts::to_symmetric_matrix( 
                flu::contacts::table_to_contacts(polymod_data,
                    age_group_limits ), 
                age_data );

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

//' Run the SEIR model for the given parameters
//'
//' @param population The population size of the different age groups, subdivided into risk groups 
//' @param initial_infected The corresponding number of initially infected
//' @param vaccine_calendar A vaccine calendar valid for that year
//' @param contact_matrix Contact rates between different age groups
//' @param susceptibility Vector with susceptibilities of each age group
//' @param transmissibility The transmissibility of the strain
//' @param infection_delays Vector with the time of latent infection and time infectious
//' @param dates Dates to return values for.
//' @return A data frame with number of new cases after each interval during the year
//'
// [[Rcpp::export(name="infectionODEs.cpp")]]
Rcpp::DataFrame infectionODEs(
        Rcpp::NumericVector population,
        Eigen::VectorXd initial_infected, 
        flu::vaccine::vaccine_t vaccine_calendar,
        Eigen::MatrixXd contact_matrix,
        Eigen::VectorXd susceptibility, 
        double transmissibility, 
        Eigen::VectorXd infection_delays, 
        Rcpp::DateVector dates )
{

    Eigen::VectorXd popv(population.size());
    for (int i = 0; i < population.size(); ++i)
        popv[i] = population[i];
    
    if (contact_matrix.cols() != contact_matrix.rows()) 
        ::Rf_error("Contact matrix should be a square matrix");
    else if (contact_matrix.cols() != susceptibility.size())
        ::Rf_error("Contact matrix and susceptibility vector should use the same number of age groups.");
    else if (contact_matrix.cols()*3 < popv.size()) 
        ::Rf_error("Maximum of three risk groups are expected. Population vector should have the initial popv of each group");
    else if (popv.size()%contact_matrix.cols()!=0)
        ::Rf_error("Population groups and contact_matrix size mismatch");
    else if (popv.size() != initial_infected.size())
        ::Rf_error("Population vector and initial_infected should have the same number of entries");

    auto dim = popv.size();
    if (contact_matrix.cols() != popv.size()/3)
    {
        popv.conservativeResize(contact_matrix.cols()*3);
        initial_infected.conservativeResize(contact_matrix.cols()*3);
        for( size_t i = dim; i<popv.size(); ++i)
        {
            popv[i] = 0;
            initial_infected[i] = 0;
        }
    }

    std::vector<boost::posix_time::ptime> datesC;
    for( auto & d : dates )
    {
        datesC.push_back( date_to_ptime( d ) );
    }

    auto result = flu::infectionODE(
        popv, initial_infected, 
        infection_delays[0], infection_delays[1],
        susceptibility, contact_matrix, transmissibility,
        vaccine_calendar, datesC );

    Rcpp::List resultList( dim + 1 );
    Rcpp::CharacterVector columnNames;

    //Rcpp::DataFrame densities = Rcpp::wrap<Rcpp::DataFrame>( result.cases );
    // Convert times
    auto times = Rcpp::DateVector(result.times.size());
    for ( size_t i = 0; i < result.times.size(); ++i )
    {
        times[i] = ptime_to_date( result.times[i] );
        if (dates[i+1]!=times[i])
        {
            ::Rf_error("Dates do not match");
        }
    }

    if (population.hasAttribute("names"))
        columnNames = population.attr("names");

    columnNames.push_front( "Time" );
    resultList[0] = times;

    for (int i=0; i<dim; ++i)
    {
        resultList[i+1] = Eigen::VectorXd(result.cases.col(i));
        if (!population.hasAttribute("names"))
            columnNames.push_back( 
                "V" + boost::lexical_cast<std::string>( i+1 ) );
    }

    auto df = Rcpp::DataFrame(resultList);

    df.attr("names") = columnNames;

    return df;
}

//' Returns log likelihood of the predicted number of cases given the data for that week
//'
//' The model results in a prediction for the given number of new cases in a certain age group and for a certain week. This function calculates the likelihood of that given the data on reported Influenza Like Illnesses and confirmed samples.
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
//' @seealso{\link{total_log_likelihood_cases}}
//'
// [[Rcpp::export(name=".log_likelihood_cases")]]
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

//' Returns log likelihood of the predicted number of cases given the data
//'
//' The model results in a prediction for the number of new cases in a certain age group and for a certain week. This function sum the log likelihood for the predicted cases for each week and age group given the data on reported Influenza Like Illnesses and confirmed samples.
//'
//' @param epsilon Parameter for the probability distribution by age group
//' @param psi Parameter for the probability distribution
//' @param predicted Number of cases predicted by your model for each week and age group
//' @param population_size The total population size in the age groups 
//' @param ili_cases The number of Influenza Like Illness cases by week and age group
//' @param ili_monitored The size of the population monitored for ILI  by week and age group
//' @param confirmed_positive The number of samples positive for the Influenza strain  by week and age group
//' @param confirmed_samples Number of samples tested for the Influenza strain  by week and age group
//'
//'
// [[Rcpp::export(name="log_likelihood_cases")]]
double total_log_likelihood(  Eigen::VectorXd epsilon, double psi, 
        Eigen::MatrixXi predicted, Eigen::VectorXi population_size, 
        Eigen::MatrixXi ili_cases, Eigen::MatrixXi ili_monitored,
        Eigen::MatrixXi confirmed_positive, Eigen::MatrixXi confirmed_samples )
{
    double ll = 0;
    for (size_t j = 0; j < ili_cases.cols(); ++j) {
        for (size_t i = 0; i < ili_cases.rows(); ++i) {
            ll += flu::log_likelihood( epsilon[j], psi,
                    predicted(i,j), population_size[j],
                    ili_cases(i,j), ili_monitored(i,j),
                    confirmed_positive(i,j), confirmed_samples(i,j), 
                    2 );
        }
    }
    return ll;
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
//' @param outfun A function that is called for each batch. Can be useful to log certain values. 
//' @param acceptfun A function that is called whenever a sample is accepted. 
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
        Rcpp::Function outfun,
        Rcpp::Function acceptfun,
        size_t nburn,
        Eigen::VectorXd initial, 
        size_t nbatch, size_t blen = 1, bool verbose = false )
{
    auto cppLprior = [&lprior]( const Eigen::VectorXd &pars ) {
        double lPrior = Rcpp::as<double>(lprior( pars ));
        return lPrior;
    };

    auto cppLlikelihood = [&llikelihood]( const Eigen::VectorXd &pars ) {
        double ll = Rcpp::as<double>(llikelihood( pars ));
        return ll;
    };

    auto mcmcResult = flu::adaptiveMCMC( cppLprior, cppLlikelihood, outfun, acceptfun,
            nburn, initial, nbatch, blen, verbose );
    Rcpp::List rState;
    rState["batch"] = Rcpp::wrap( mcmcResult.batch );
    rState["llikelihoods"] = Rcpp::wrap( mcmcResult.llikelihoods );
    return rState;
}

//' Create a contact matrix based on polymod data.
//'
//' @param polymod_data Contact data for different age groups
//' @param demography A vector with the population size by each age {0,1,2,..}
//' @param age_group_limits The upper limits of the different age groups (by default: c(1,5,15,25,45,65), which corresponds to age groups: <1, 1-14, 15-24, 25-44, 45-64, >=65.
//'
//' @return Returns a symmetric matrix with the frequency of contact between each age group
//'
// [[Rcpp::export]]
Eigen::MatrixXd contact_matrix(  
        Eigen::MatrixXi polymod_data,
        std::vector<size_t> demography,
        Rcpp::NumericVector age_group_limits = Rcpp::NumericVector::create(
            1, 5, 15, 25, 45, 65 ) )
{

    if (polymod_data.cols() - 2 != age_group_limits.size() + 1)
        ::Rf_error("Number of age groups should be consistent for the polymod_data and the age_group_limits");
    
    flu::data::age_data_t age_data;
    age_data.age_sizes = demography;

    auto agl_v = std::vector<size_t>( 
                age_group_limits.begin(), age_group_limits.end() );
    age_data.age_group_sizes = flu::data::group_age_data( demography, agl_v );

    return flu::contacts::to_symmetric_matrix( 
            flu::contacts::table_to_contacts(polymod_data, agl_v), 
            age_data );
}

size_t as_age_group_index( size_t age,
        Rcpp::NumericVector limits = Rcpp::NumericVector::create(
            1, 5, 15, 25, 45, 65 ) )
{
    std::deque<size_t> group_barriers = std::deque<size_t>( 
                limits.begin(), limits.end() );
    size_t ag = 1;
    while (group_barriers.size() != 0 
            && age >= group_barriers[0] )
    {
        ++ag;
        group_barriers.pop_front();
    }
    return ag;
}

//' Create age group level description based on passed upper limits
//'
//' @description Returns a vector of age group levels given the upper age group limits. These levels can be used as the named levels in a factor
//'
//' @param limits The upper limit to each age groups (not included) (1,5,15,25,45,65) corresponds to the following age groups: <1, 1-4, 5-14, 15-24, 25-44, 45-64 and >=65.
//'
//' @return Vector representing the age group(s)
//'
// [[Rcpp::export]]
Rcpp::CharacterVector age_group_levels(Rcpp::NumericVector limits = Rcpp::NumericVector::create()) {
    Rcpp::CharacterVector lvls;
    // Create levels
    for (auto i = 0; i <= limits.size(); ++i) {
        if (i == 0) {
            Rcpp::String str("[0,");
            if (limits.size() != 0) {
                str += limits[i];
            } else {
                str += "+";
            }
            str += ")";
            lvls.push_back(str);
        } else if (i == limits.size()) {
            Rcpp::String str("[");
            str += limits[i-1];
            str += ",+)";
            lvls.push_back(str);
        } else {
            Rcpp::String str("[");
            str += limits[i-1];
            str += ",";
            str += limits[i];
            str += ")";
            lvls.push_back(str);
        }
    }

    return lvls;
}

//' Extract upper age group limits from age group level description
//'
//' @description Returns a vector of age group limits given the age group level descriptions. This is a helper function, which is essentially the reverse of \code{\link{age_group_levels}}.
//'
//' @param levels The levels representing each age groups.
//'
//' @return Vector representing the age group(s) limits
//'
//' @seealso \code{\link{age_group_levels}} For the reverse of this function.
//'
// [[Rcpp::export]]
Rcpp::IntegerVector age_group_limits(std::vector<std::string> levels) 
{
    //Rcpp::IntegerVector age_group_limits(Rcpp::CharacterVector levels) 
    Rcpp::IntegerVector limits;

    // Regex support is broken in gcc 4.8, so we will scan by hand
    /*auto r = std::regex("\\d+,\\s*(\\d+)");
    std::smatch match;
    for(auto&& lvl : levels) {
        if (regex_search(lvl, match, r)) {
            limits.push_back(std::stoi(match[1]));
        }
    }
    */

    for(auto&& lvl : levels) {
        if (lvl.length() > 1) {
            std::string match;
            for(int i = (int)lvl.length()-1; i >= 0; --i) {
                if (isdigit(lvl[i])) {
                    match = lvl[i] + match;
                } if (lvl[i] == ',') {
                    if(match.length() > 0)
                        limits.push_back(std::stoi(match));
                    break;
                }
            }
        }
    }

    return limits;
}

//' Age as age group
//'
//' @description Returns the age group a certain age belongs to given the upper age group limits 
//'
//' @param age The relevant age. This can be a vector.
//' @param limits The upper limit to each age groups (not included) (1,5,15,25,45,65) corresponds to the following age groups: <1, 1-4, 5-14, 15-24, 25-44, 45-64 and >=65.
//'
//' @return Factors representing the age group(s)
//'
// [[Rcpp::export]]
Rcpp::IntegerVector as_age_group( Rcpp::NumericVector age,
        Rcpp::NumericVector limits = Rcpp::NumericVector::create(
            1, 5, 15, 25, 45, 65 ) )
{
    Rcpp::IntegerVector out;
    Rcpp::CharacterVector lvls = age_group_levels(limits);

    // Map to indexes
    for (auto &&ag : age)
        out.push_back(as_age_group_index(ag, limits));
    out.attr("class") = "factor";
    out.attr("levels") = lvls;
    return out;
}

//' @title Stratify the population by age
//'
//' @description Stratifies the population and returns the population size of each age group.
//'
//' @param age_sizes A vector containing the population size by age (first element is number of people of age 1 and below)
//' @param limits The upper limit to each age groups (not included) (1,5,15,25,45,65) corresponds to the following age groups: <1, 1-4, 5-14, 15-24, 25-44, 45-64 and >=65.
//'
//' @return A vector with the population in each age group.
//'
// [[Rcpp::export(name="stratify_by_age")]]
Rcpp::IntegerVector separate_into_age_groups( std::vector<size_t> age_sizes,
        Rcpp::NumericVector limits = Rcpp::NumericVector::create(
            1, 5, 15, 25, 45, 65 ) )
{
    Rcpp::CharacterVector lvls = age_group_levels(limits);
    auto agl_v = std::vector<size_t>( 
                limits.begin(), limits.end() );

    auto v = flu::data::group_age_data( age_sizes, agl_v );
    Rcpp::IntegerVector out = Rcpp::IntegerVector(v.data(), v.data() + v.size());
    out.attr("names") = lvls;
    return out;
}

//' @title Stratify age groups into different risk groups
//' 
//' @description Stratifies the age groups and returns the population size of each age group and risk group.
//'
//' @param age_groups A vector containing the population size of each age group
//' @param risk A matrix with the fraction in the risk groups. The leftover fraction is assumed to be low risk
//'
//' @return A vector with the population in the low risk groups, followed by the other risk groups. The length is equal to the number of age groups times the number of risk groups (including the low risk group).
//'
// [[Rcpp::export(name=".stratify_by_risk")]]
Eigen::VectorXd stratify_by_risk( 
        const Eigen::VectorXd &age_groups, const Eigen::VectorXd &risk, size_t no_risk_groups )
{
    return flu::data::stratify_by_risk(age_groups, risk, no_risk_groups);
}

//' @title Calculate R0 from transmission rate
//'
//' @description Uses the transmission rate (\eqn{\lambda}), contact matrix (\eqn{c}), population (\eqn{N}), and infectious period (\eqn{\gamma}) 
//' to calculate the R0 using the following equation.
//' \deqn{\lambda max(EV(C)) \gamma}
//' where \eqn{EV(C)} denotes the eigenvalues of the matrix \eqn{C} which is calculated from the contact matrix and the population 
//' (\eqn{C[i,j] = c[i,j] N[j]}).
//'
//' @param transmission_rate The transmission rate of the disease
//' @param contaxt_matrix The contact matrix between age groups
//' @param age_groups The population size of the different age groups
//' @param duration Duration of the infectious period. Default value is 1.8 days
//'
//' @return Returns the R0
// [[Rcpp::export]]
double as_R0(double transmission_rate, Eigen::MatrixXd contact_matrix, Eigen::VectorXd age_groups,
        double duration = 1.8) {
    //auto evs = (contact_matrix*population).eigenvalues();
    //return transmission_rate*evs.maxCoeff()*duration;
    if (contact_matrix.cols() != age_groups.size())
        ::Rf_error("Number of age groups should be equal to the dimensions of the contact_matrix");
    auto a = contact_matrix;
    for (size_t i = 0; i < a.rows(); ++i) {
        for (size_t j = 0; j < a.cols(); ++j) {
            a(i,j) = contact_matrix(i,j)*age_groups[j];
        }
    }

    auto evs = a.eigenvalues().real().maxCoeff();
    return transmission_rate*evs*duration;
}

//' @title Calculate transmission rate from R0 
//'
//' @description Uses the R0 (\eqn{R0}), contact matrix (\eqn{c}), population (\eqn{N}), and infectious period (\eqn{\gamma}) 
//' to calculate the transmission rate using the following equation.
//' \deqn{R0/(max(EV(C)) \gamma)}
//' where \eqn{EV(C)} denotes the eigenvalues of the matrix \eqn{C} which is calculated from the contact matrix and the population 
//' (\eqn{C[i,j] = c[i,j] N[j]}).
//'
//' @param R0 The R0 of the disease
//' @param contaxt_matrix The contact matrix between age groups
//' @param age_groups The population size of the different age groups
//' @param duration Duration of the infectious period. Default value is 1.8 days
//'
//' @return Returns the transmission rate 
// [[Rcpp::export]]
double as_transmission_rate(double R0, Eigen::MatrixXd contact_matrix, Eigen::VectorXd age_groups,
        double duration = 1.8) {
    //auto evs = (contact_matrix*population).eigenvalues();
    //return transmission_rate*evs.maxCoeff()*duration;
    if (contact_matrix.cols() != age_groups.size())
        ::Rf_error("Number of age groups should be equal to the dimensions of the contact_matrix");
    auto a = contact_matrix;
    for (size_t i = 0; i < a.rows(); ++i) {
        for (size_t j = 0; j < a.cols(); ++j) {
            a(i,j) = contact_matrix(i,j)*age_groups[j];
        }
    }

    auto evs = a.eigenvalues().real().maxCoeff();
    return R0/(evs*duration);
}
