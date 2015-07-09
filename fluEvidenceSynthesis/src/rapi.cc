#include <boost/date_time.hpp>

#include "rcppwrap.hh"

#include "proposal.hh"
#include "contacts.hh"
#include "model.hh"

namespace bt = boost::posix_time;

/**
 * This file contains some thin wrappers around C++ code. Mostly used
 * for testing the C++ code from R.
 *
 * The wrappers are needed, because in the C++ code we often want to pass
 * by const reference, which is impossible in a function called from R
 */

// [[Rcpp::export]]
Eigen::VectorXd updateMeans( Eigen::VectorXd means,
        Eigen::VectorXd v, size_t n )
{
    return flu::proposal::updateMeans( means, v, n );
}

// [[Rcpp::export]]
Eigen::MatrixXd updateCovariance( Eigen::MatrixXd cov, 
        Eigen::VectorXd v, Eigen::VectorXd means, size_t n )
{
    return flu::proposal::updateCovariance( cov, v, means, n );
}

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

// [[Rcpp::export]]
Rcpp::DataFrame runSEIRModel(
        std::vector<size_t> age_sizes, 
        flu::vaccine::vaccine_t vaccine_calendar,
        flu::contacts::contacts_t polymod_data,
        flu::state_t current_state )
{

    double result[7644];

    flu::data::age_data_t age_data;
    age_data.age_sizes = age_sizes;
    age_data.age_group_sizes = flu::data::group_age_data( age_sizes );
    auto pop_vec = flu::data::separate_into_risk_groups( 
            age_data.age_group_sizes );

    /*translate into an initial infected population*/
    double curr_init_inf[NAG];
    for(int i=0;i<NAG;i++)
        curr_init_inf[i]=pow(10,current_state.parameters.init_pop);

    auto curr_c = flu::contacts::shuffle_by_id( polymod_data, 
            current_state.contact_ids );

    auto current_contact_regular = 
        flu::contacts::to_symmetric_matrix( curr_c, age_data );

    flu::one_year_SEIR_with_vaccination(result, pop_vec, curr_init_inf, current_state.time_latent, current_state.time_infectious, current_state.parameters.susceptibility, current_contact_regular, current_state.parameters.transmissibility, vaccine_calendar );

    //if((((int)(t*step_rate))%step_rate)==step_rate/2)
    auto resultMatrix = Eigen::MatrixXd( 7644/21, 21 );
    for(size_t t = 0; t<resultMatrix.rows(); ++t)
    {
        for(size_t i=0;i<NAG;i++)
        {
            resultMatrix(t, i) = result[t*3*NAG+i];
            resultMatrix(t, i+NAG) = result[t*3*NAG+NAG+i];
            resultMatrix(t, i+2*NAG) = result[t*3*NAG+2*NAG+i];
        }
    }
    Rcpp::DataFrame densities = Rcpp::wrap<Rcpp::DataFrame>( resultMatrix );
    return densities;
    //return resultMatrix;
}
