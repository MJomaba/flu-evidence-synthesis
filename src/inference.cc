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
std::vector<state_t> inference( std::vector<size_t> age_sizes, 
        Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
        Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
        flu::vaccine::vaccine_t vaccine_calendar,
        flu::contacts::contacts_t polymod_data,
        flu::state_t init_sample, 
        int mcmc_chain_length = 100000, 
        int burn_in = 10000, int thinning = 100 )
{
    std::vector<state_t> results;

    double prop_init_inf[NAG];
    double curr_init_inf[NAG];
    int step_mat;
    double pop_RCGP[5];
    double prop_likelihood;

    double my_acceptance_rate;

    state_t current_state = init_sample;

    flu::data::age_data_t age_data;
    age_data.age_sizes = age_sizes;
    age_data.age_group_sizes = flu::data::group_age_data( age_sizes );
    auto pop_vec = flu::data::separate_into_risk_groups( 
            age_data.age_group_sizes );

    /*pop RCGP*/
    pop_RCGP[0]=pop_vec[0]+pop_vec[1]+pop_vec[7]+pop_vec[8]+pop_vec[14]+pop_vec[15];
    pop_RCGP[1]=pop_vec[2]+pop_vec[9]+pop_vec[16];
    pop_RCGP[2]=pop_vec[3]+pop_vec[4]+pop_vec[10]+pop_vec[11]+pop_vec[17]+pop_vec[18];
    pop_RCGP[3]=pop_vec[5]+pop_vec[12]+pop_vec[19];
    pop_RCGP[4]=pop_vec[6]+pop_vec[13]+pop_vec[20];

    /********************************************************************************************************************************************************
    initialisation point to start the MCMC
    *********************************************************************************************************************************************************/

    /*translate into an initial infected population*/
    for(int i=0;i<NAG;i++)
        curr_init_inf[i]=pow(10,current_state.parameters.init_pop);


    auto curr_c = contacts::shuffle_by_id( polymod_data, 
            current_state.contact_ids );

    auto current_contact_regular = 
        contacts::to_symmetric_matrix( curr_c, age_data );

    auto result = one_year_SEIR_with_vaccination(pop_vec, curr_init_inf, current_state.time_latent, current_state.time_infectious, current_state.parameters.susceptibility, current_contact_regular, current_state.parameters.transmissibility, vaccine_calendar, 7*24 );

    /*curr_psi=0.00001;*/
    auto d_app = 3;
    current_state.likelihood = log_likelihood_hyper_poisson(
            current_state.parameters.epsilon, current_state.parameters.psi, 
            days_to_weeks_5AG( result ), 
            ili, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

    auto proposal_state = proposal::initialize( 9 );

    /**************************************************************************************************************************************************************
    Start of the MCMC
    **************************************************************************************************************************************************************/

    step_mat=1;            /*number of contacts exchanged*/
    double p_ac_mat=0.10;          /*prob to redraw matrices*/

    for(int k=1; k<=mcmc_chain_length + burn_in; k++)
    {

        /*update of the variance-covariance matrix and the mean vector*/
        proposal_state = proposal::update( std::move( proposal_state ),
                current_state.parameters, k );

        if(k%thinning==0 && k>burn_in)
        {
            // Add results
            for( size_t i = 0; i < POLY_PART; ++i )
                current_state.contact_ids[i] = curr_c.contacts[i].id;

            results.push_back(current_state);
        }

        auto proposed_par = proposal::haario_adapt_scale(
                current_state.parameters,
                proposal_state.chol_emp_cov,
                proposal_state.chol_ini,100,0.05, 
                proposal_state.adaptive_scaling );

        auto prior_ratio = 
            log_prior(proposed_par, current_state.parameters, false );

        if (!std::isfinite(prior_ratio))
        {
            //Rcpp::Rcout << "Invalid proposed par" << std::endl;
            // TODO: code duplication with failure of acceptance
            if(k>=1000)
                proposal_state.adaptive_scaling
                    -=0.234*proposal_state.conv_scaling;
        } else {
            /*translate into an initial infected population*/
            for(int i=0;i<NAG;i++)
                prop_init_inf[i]=pow(10,proposed_par.init_pop);

            auto prop_c = curr_c;

            /*do swap of contacts step_mat times (reduce or increase to change 'distance' of new matrix from current)*/
            // TODO/WARN Need to draw this before hand and pass it as data to
            // likelihood function... Even when doing that we still need to know k,
            // so might as well make the likelihood function increase k when called
            if(R::runif(0,1) < p_ac_mat)
                prop_c = contacts::bootstrap_contacts( std::move(prop_c),
                        polymod_data, step_mat );

            auto prop_contact_regular = 
                contacts::to_symmetric_matrix( prop_c, age_data );

            result = one_year_SEIR_with_vaccination(pop_vec, prop_init_inf, current_state.time_latent, current_state.time_infectious, proposed_par.susceptibility, prop_contact_regular, proposed_par.transmissibility, vaccine_calendar, 7*24 );

            /*computes the associated likelihood with the proposed values*/
            prop_likelihood=log_likelihood_hyper_poisson(
                    proposed_par.epsilon, proposed_par.psi, 
                    days_to_weeks_5AG(result), ili, mon_pop, 
                    n_pos, n_samples, pop_RCGP, d_app);

            /*Acceptance rate include the likelihood and the prior but no correction for the proposal as we use a symmetrical RW*/
            // Make sure accept works with -inf prior
            // MCMC-R alternative prior?
            if (std::isinf(prop_likelihood) && std::isinf(current_state.likelihood) )
                my_acceptance_rate = exp(prior_ratio); // We want to explore and find a non infinite likelihood
            else 
                my_acceptance_rate=
                    exp(prop_likelihood-current_state.likelihood+
                    prior_ratio);

            if(R::runif(0,1)<my_acceptance_rate) /*with prior*/
            {
                /*update the acceptance rate*/
                proposal_state.acceptance++;
                if(k>=1000)
                    proposal_state.adaptive_scaling
                        += 0.766*proposal_state.conv_scaling;

                current_state.parameters = proposed_par;

                /*update current likelihood*/
                current_state.likelihood=prop_likelihood;

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
    }
    return results;
}

namespace flu {
    mcmc_result_t binomial_inference( std::vector<size_t> age_sizes, 
            Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
            Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
            flu::vaccine::vaccine_t vaccine_calendar,
            flu::contacts::contacts_t polymod_data,
            Eigen::VectorXd initial, 
            int nburn = 10000, 
            int nbatch = 1000, int blen = 1 )
    {
        mcmc_result_t results;
        results.batch = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>( nbatch, initial.size() );
        results.llikelihoods = Eigen::VectorXd( nbatch );


        double prop_init_inf[NAG];
        double curr_init_inf[NAG];
        double pop_RCGP[5];

        double my_acceptance_rate;

        flu::data::age_data_t age_data;
        age_data.age_sizes = age_sizes;
        age_data.age_group_sizes = flu::data::group_age_data( age_sizes );
        auto pop_vec = flu::data::separate_into_risk_groups( 
                age_data.age_group_sizes );

        /*pop RCGP*/
        pop_RCGP[0]=pop_vec[0]+pop_vec[1]+pop_vec[7]+pop_vec[8]+pop_vec[14]+pop_vec[15];
        pop_RCGP[1]=pop_vec[2]+pop_vec[9]+pop_vec[16];
        pop_RCGP[2]=pop_vec[3]+pop_vec[4]+pop_vec[10]+pop_vec[11]+pop_vec[17]+pop_vec[18];
        pop_RCGP[3]=pop_vec[5]+pop_vec[12]+pop_vec[19];
        pop_RCGP[4]=pop_vec[6]+pop_vec[13]+pop_vec[20];

        /********************************************************************************************************************************************************
          initialisation point to start the MCMC
         *********************************************************************************************************************************************************/

        auto current_pars = initial;
        /*translate into an initial infected population*/
        for(int i=0;i<NAG;i++)
            curr_init_inf[i]=pow(10,current_pars[7] );

        std::vector<size_t> current_contact_ids;
        current_contact_ids.reserve( polymod_data.contacts.size() );
        for ( size_t i = 0; i < polymod_data.contacts.size(); ++i )
            current_contact_ids.push_back( i );


        auto time_latent = 0.8;
        auto time_infectious = 1.8;
        auto curr_c = contacts::shuffle_by_id( polymod_data, 
                current_contact_ids );

        auto current_contact_regular = 
            contacts::to_symmetric_matrix( curr_c, age_data );

        Eigen::VectorXd susc(3);
        susc << current_pars[4], current_pars[5], current_pars[6];

        auto result = one_year_SEIR_with_vaccination(pop_vec, 
                curr_init_inf, time_latent, time_infectious, 
                susc,
                current_contact_regular, current_pars[3], 
                vaccine_calendar, 7*24 );

        /*curr_psi=0.00001;*/
        long double current_llikelihood = binomial_log_likelihood_year(
                { current_pars[0], current_pars[0], 
                current_pars[1], current_pars[1], current_pars[2] },
                days_to_weeks_5AG( result ), 
                ili, mon_pop, n_pos, n_samples, pop_RCGP);

        auto log_prior = []( Eigen::VectorXd pars )
        {
            if (pars[0] <= 0 || pars[0] >= 1 ||
                    pars[1] <= 0 || pars[1] >= 1 ||
                    pars[2] <= 0 || pars[2] >= 1 ||
                    pars[3] <= 0 ||
                    pars[4] <= 0 || pars[4] >= 1 ||
                    pars[5] <= 0 || pars[5] >= 1 ||
                    pars[6] <= 0 || pars[6] >= 1 ||
                    pars[7] < log(0.00001) || pars[7] > log(10))
                return log(0);

            return R::dnorm(pars[3], 0.1653183, 0.02773, 1);
        };

        auto current_lprior = log_prior( current_pars );

        auto proposal_state = proposal::initialize( current_pars.size() );

        /**************************************************************************************************************************************************************
          Start of the MCMC
         **************************************************************************************************************************************************************/

        int step_mat=1;            /*number of contacts exchanged*/
        double p_ac_mat=0.10;          /*prob to redraw matrices*/

        size_t sampleCount = 0;
        int k = 0;
        while( sampleCount < nbatch )
        {
            ++k;

            /*update of the variance-covariance matrix and the mean vector*/
            proposal_state = proposal::update( std::move( proposal_state ),
                    current_pars, k );

            auto proposed_pars = proposal::haario_adapt_scale(
                    current_pars,
                    proposal_state.chol_emp_cov,
                    proposal_state.chol_ini,100,0.05, 
                    proposal_state.adaptive_scaling );

            auto proposed_lprior = log_prior( proposed_pars );

            /*translate into an initial infected population*/
            for(int i=0;i<NAG;i++)
                prop_init_inf[i]=pow(10,proposed_pars[7]);

            auto prop_c = curr_c;

            if(R::runif(0,1) < p_ac_mat)
                prop_c = contacts::bootstrap_contacts( std::move(prop_c),
                        polymod_data, step_mat );

            auto prop_contact_regular = 
                contacts::to_symmetric_matrix( prop_c, age_data );

            susc << proposed_pars[4], proposed_pars[5], proposed_pars[6];
            result = one_year_SEIR_with_vaccination(pop_vec, 
                    curr_init_inf, time_latent, time_infectious, 
                    susc,
                    prop_contact_regular, proposed_pars[3], 
                    vaccine_calendar, 7*24 );

            long double proposed_llikelihood = binomial_log_likelihood_year(
                    { proposed_pars[0], proposed_pars[0], 
                    proposed_pars[1], proposed_pars[1], proposed_pars[2] },
                    days_to_weeks_5AG( result ), 
                    ili, mon_pop, n_pos, n_samples, pop_RCGP);

            /*Acceptance rate include the likelihood and the prior but no correction for the proposal as we use a symmetrical RW*/
            // Make sure accept works with -inf prior
            // MCMC-R alternative prior?
            if (std::isinf(proposed_llikelihood) 
                    && std::isinf(current_llikelihood) )
                my_acceptance_rate = exp(proposed_lprior-current_lprior); // We want to explore and find a non infinite likelihood
            else 
                my_acceptance_rate=
                    exp(proposed_llikelihood - current_llikelihood +
                            proposed_lprior - current_lprior);

            if(R::runif(0,1)<my_acceptance_rate) /*with prior*/
            {
                /*update the acceptance rate*/
                proposal_state.acceptance++;
                if(k>=1000)
                    proposal_state.adaptive_scaling
                        += 0.766*proposal_state.conv_scaling;

                current_pars = proposed_pars;

                /*update current likelihood*/
                current_llikelihood = proposed_llikelihood;
                current_lprior = proposed_lprior;

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
            if(k%blen==0 && k>=nburn)
            {
                // Add results
                results.llikelihoods[sampleCount] = current_llikelihood;
                results.batch.row( sampleCount ) = current_pars;
                ++sampleCount;
            }
        }
        return results;
    }
}
