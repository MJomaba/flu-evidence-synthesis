#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <iostream>

#include "model.hh"
#include "state.hh"
#include "data.hh"
#include "contacts.hh"
#include "vaccine.hh"
#include "proposal.hh"

#include "rcppwrap.hh"
#include<RcppEigen.h>

using namespace flu;

// [[Rcpp::export]]
std::vector<state_t> inference( std::vector<size_t> age_sizes, 
        Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
        Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
        flu::vaccine::vaccine_t vaccine_calendar,
        flu::contacts::contacts_t polymod_data,
        flu::state_t init_state, 
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

    state_t current_state = init_state;

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

    auto result = one_year_SEIR_with_vaccination(pop_vec, curr_init_inf, current_state.time_latent, current_state.time_infectious, current_state.parameters.susceptibility, current_contact_regular, current_state.parameters.transmissibility, vaccine_calendar );

    /*curr_psi=0.00001;*/
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
            if(R::runif(0,1) < p_ac_mat)
                prop_c = contacts::bootstrap_contacts( std::move(prop_c),
                        polymod_data, step_mat );

            auto prop_contact_regular = 
                contacts::to_symmetric_matrix( prop_c, age_data );

            result = one_year_SEIR_with_vaccination(pop_vec, prop_init_inf, current_state.time_latent, current_state.time_infectious, proposed_par.susceptibility, prop_contact_regular, proposed_par.transmissibility, vaccine_calendar );

            /*computes the associated likelihood with the proposed values*/
            prop_likelihood=log_likelihood_hyper_poisson(
                    proposed_par.epsilon, proposed_par.psi, 
                    days_to_weeks_5AG(result), ili, mon_pop, 
                    n_pos, n_samples, pop_RCGP, d_app);

            /*Acceptance rate include the likelihood and the prior but no correction for the proposal as we use a symmetrical RW*/
            // Make sure accept works with -inf prior
            // MCMC-R alternative prior?
            my_acceptance_rate=exp(prop_likelihood-current_state.likelihood+
                    log_prior(proposed_par, current_state.parameters, false ));

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
