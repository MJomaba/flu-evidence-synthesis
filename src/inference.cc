#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <iostream>

#include "model.h"
#include "state.h"
#include "data.h"
#include "contacts.h"
#include "vaccine.h"
#include "proposal.h"

#include "mcmc.h"

#include "rcppwrap.h"
#include<RcppEigen.h>

using namespace flu;

//' MCMC based inference of the parameter values given the different data sets
//'
//' @param age_sizes A vector with the population size by each age {1,2,..}
//' @param ili The number of Influenza-like illness cases per week
//' @param mon_pop The number of people monitored for ili
//' @param n_pos The number of positive samples for the given strain (per week)
//' @param n_samples The total number of samples tested 
//' @param vaccine_calendar A vaccine calendar valid for that year
//' @param polymod_data Contact data for different age groups
//' @param init_sample The initial parameters needed to run the ODE model (typically one of the posterior sample created when running the inference)
//' @param mcmc_chain_length The number of MCMC steps to sample from
//' @param thinning Keep every so many samples
//' @param burn_in The number of initial samples to skip
//' @return A vector with posterior samples of the parameters (length of \code{mcmc_chain_length}/\code{thinning})
//'
//'
// [[Rcpp::export]]
mcmc_result_inference_t inference( std::vector<size_t> age_sizes, 
        Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
        Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
        flu::vaccine::vaccine_t vaccine_calendar,
        Eigen::MatrixXi polymod_data,
        Eigen::VectorXd initial, 
        size_t nburn = 0,
        size_t nbatch = 1000, size_t blen = 1 )
{
    mcmc_result_inference_t results;
    results.batch = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>( nbatch, initial.size() );
    results.llikelihoods = Eigen::VectorXd( nbatch );
    results.contact_ids = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>( nbatch, polymod_data.rows() );


    auto nag = 7;

    int step_mat;
    double pop_RCGP[5];
    double prop_likelihood;

    double my_acceptance_rate;

    std::vector<size_t> age_group_limits = {1,5,15,25,45,65};

    flu::data::age_data_t age_data;
    age_data.age_sizes = age_sizes;
    age_data.age_group_sizes = flu::data::group_age_data( age_sizes,
            age_group_limits );

    Eigen::MatrixXd risk_proportions = Eigen::MatrixXd( 
            2, age_data.age_group_sizes.size() );
    risk_proportions << 
        0.021, 0.055, 0.098, 0.087, 0.092, 0.183, 0.45, 
        0, 0, 0, 0, 0, 0, 0;

    auto pop_vec = flu::data::separate_into_risk_groups( 
            age_data.age_group_sizes, risk_proportions  );

    /*pop RCGP*/
    pop_RCGP[0]=pop_vec[0]+pop_vec[1]+pop_vec[7]+pop_vec[8]+pop_vec[14]+pop_vec[15];
    pop_RCGP[1]=pop_vec[2]+pop_vec[9]+pop_vec[16];
    pop_RCGP[2]=pop_vec[3]+pop_vec[4]+pop_vec[10]+pop_vec[11]+pop_vec[17]+pop_vec[18];
    pop_RCGP[3]=pop_vec[5]+pop_vec[12]+pop_vec[19];
    pop_RCGP[4]=pop_vec[6]+pop_vec[13]+pop_vec[20];

    /********************************************************************************************************************************************************
    initialisation point to start the MCMC
    *********************************************************************************************************************************************************/

    auto pars_to_epsilon = []( const Eigen::VectorXd &pars )
    {
        Eigen::VectorXd eps(5);
        eps << pars[0], pars[0], pars[1], pars[1], pars[2];
        return eps;
    };

    auto pars_to_susceptibility = []( const Eigen::VectorXd &pars )
    {
        Eigen::VectorXd susc(7);
        susc << pars[5], pars[5], pars[5], 
             pars[6], pars[6], pars[6], pars[7];
        return susc;
    };

    /*translate into an initial infected population*/
    std::vector<size_t> contact_ids;
    for (size_t i = 0; i < (size_t)polymod_data.rows(); ++i)
        contact_ids.push_back(i);

    auto curr_parameters = initial;
    auto curr_init_inf = Eigen::VectorXd::Constant( 
            nag, pow(10,curr_parameters[8]) );

    auto polymod = flu::contacts::table_to_contacts( polymod_data, 
            age_group_limits ); 

    auto curr_c = contacts::shuffle_by_id( polymod, 
            contact_ids );

    auto current_contact_regular = 
        contacts::to_symmetric_matrix( curr_c, age_data );

    auto time_latent = 0.8;
    auto time_infectious = 1.8;

    auto result = one_year_SEIR_with_vaccination(pop_vec, 
            curr_init_inf, 
            time_latent, time_infectious, 
            pars_to_susceptibility(curr_parameters),
            current_contact_regular, curr_parameters[4], 
            vaccine_calendar, 7*24 );

    /*curr_psi=0.00001;*/
    auto d_app = 3;
    auto curr_llikelihood = log_likelihood_hyper_poisson(
            pars_to_epsilon(curr_parameters),
            curr_parameters[3], 
            days_to_weeks_5AG( result ), 
            ili, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

    auto proposal_state = proposal::initialize( 9 );

    /**************************************************************************************************************************************************************
    Start of the MCMC
    **************************************************************************************************************************************************************/

    step_mat=1;            /*number of contacts exchanged*/
    double p_ac_mat=0.10;          /*prob to redraw matrices*/

    size_t sampleCount = 0;
    int k = 0;
    while(sampleCount<nbatch)
    {
        ++k;

        /*update of the variance-covariance matrix and the mean vector*/
        proposal_state = proposal::update( std::move( proposal_state ),
                curr_parameters, k );

        auto prop_parameters = proposal::haario_adapt_scale(
                curr_parameters,
                proposal_state.chol_emp_cov,
                proposal_state.chol_ini,100,0.05, 
                proposal_state.adaptive_scaling );

        auto prior_ratio = 
            log_prior(prop_parameters, curr_parameters, false );

        if (!std::isfinite(prior_ratio))
        {
            //Rcpp::Rcout << "Invalid proposed par" << std::endl;
            // TODO: code duplication with failure of acceptance
            if(k>=1000)
                proposal_state.adaptive_scaling
                    -=0.234*proposal_state.conv_scaling;
        } else {
            /*translate into an initial infected population*/
            
            auto prop_init_inf = Eigen::VectorXd::Constant( 
                    nag, pow(10,prop_parameters[8]) );

            auto prop_c = curr_c;

            /*do swap of contacts step_mat times (reduce or increase to change 'distance' of new matrix from current)*/
            // TODO/WARN Need to draw this before hand and pass it as data to
            // likelihood function... Even when doing that we still need to know k,
            // so might as well make the likelihood function increase k when called
            if(R::runif(0,1) < p_ac_mat)
                prop_c = contacts::bootstrap_contacts( std::move(prop_c),
                        polymod, step_mat );

            auto prop_contact_regular = 
                contacts::to_symmetric_matrix( prop_c, age_data );

            result = one_year_SEIR_with_vaccination(pop_vec, prop_init_inf, time_latent, time_infectious, 
                    pars_to_susceptibility(prop_parameters), 
                    prop_contact_regular, 
                    prop_parameters[4], vaccine_calendar, 7*24 );

            /*computes the associated likelihood with the proposed values*/
            prop_likelihood=log_likelihood_hyper_poisson(
                    pars_to_epsilon(prop_parameters), 
                    prop_parameters[3], 
                    days_to_weeks_5AG(result), ili, mon_pop, 
                    n_pos, n_samples, pop_RCGP, d_app);

            /*Acceptance rate include the likelihood and the prior but no correction for the proposal as we use a symmetrical RW*/
            // Make sure accept works with -inf prior
            // MCMC-R alternative prior?
            if (std::isinf(prop_likelihood) && std::isinf(curr_llikelihood) )
                my_acceptance_rate = exp(prior_ratio); // We want to explore and find a non infinite likelihood
            else 
                my_acceptance_rate=
                    exp(prop_likelihood-curr_llikelihood+
                    prior_ratio);

            if(R::runif(0,1)<my_acceptance_rate) /*with prior*/
            {
                /*update the acceptance rate*/
                proposal_state.acceptance++;
                if(k>=1000)
                    proposal_state.adaptive_scaling
                        += 0.766*proposal_state.conv_scaling;

                curr_parameters = prop_parameters;

                /*update current likelihood*/
                curr_llikelihood=prop_likelihood;

                /*new proposed contact matrix*/
                /*update*/
                curr_c = prop_c;
                current_contact_regular=prop_contact_regular;
            }
            else /*if reject*/
            {
                if(k>=1000)
                    proposal_state.adaptive_scaling
                        -=0.234*proposal_state.conv_scaling;
            }
        }

        if(k%blen==0 && k>=(int)nburn)
        {
            // Add results
            results.llikelihoods[sampleCount] = curr_llikelihood;
            results.batch.row( sampleCount ) = curr_parameters;
            for( size_t i = 0; i < curr_c.contacts.size(); ++i )
                results.contact_ids( sampleCount, i ) =
                    curr_c.contacts[i].id;

            ++sampleCount;
        }
    }
    return results;
}
