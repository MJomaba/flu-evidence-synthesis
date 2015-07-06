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
    double result[7644], result_by_week[260]; /*21*52 number of new cases per week*/
    double lv, prop_likelihood;

    double my_acceptance_rate;
    /*FILE *log_file;
    FILE *f_pos_sample, *f_n_sample, *f_GP, *f_mon_pop;
    FILE *f_posterior;*/

    state_t current_state = init_state;

    /*int mcmc_chain_length=10000000;
    int burn_in=1000000;
    int thinning=1000;*/

    /*opens the file with the number of positive samples for that strain and season*/
    //f_pos_sample=read_file(data_path,"positivity.txt");

    /*opens the file with the total number of samples for that strain and season*/
    //f_n_sample=read_file(data_path,"n_samples.txt");

    /*opens the file with the RCGP ILI numbers*/
    //f_GP=read_file(data_path,"ILI.txt");

    /*opens the file with the size of the monitored population*/
    //f_mon_pop=read_file(data_path,"mon_pop.txt");

    /*********************************************************************************************************************************************************************************
    Initialisation bit LOADING DATA
    **********************************************************************************************************************************************************************************/

    /*load the positivity data*/
    /*for(int i=0;i<52;i++)
        save_fscanf(f_pos_sample,"%d %d %d %d %d",&n_pos[i*5],&n_pos[i*5+1],&n_pos[i*5+2],&n_pos[i*5+3],&n_pos[i*5+4]);*/

    /*load the number of samples data*/
    /*for(int i=0;i<52;i++)
        save_fscanf(f_n_sample,"%d %d %d %d %d",&n_samples[i*5],&n_samples[i*5+1],&n_samples[i*5+2],&n_samples[i*5+3],&n_samples[i*5+4]);*/

    /*load the number of GP ILI consultation data*/
    /*for(int i=0;i<52;i++)
        save_fscanf(f_GP,"%d %d %d %d %d",&ILI[i*5],&ILI[i*5+1],&ILI[i*5+2],&ILI[i*5+3],&ILI[i*5+4]);*/

    /*load the size of the monitored population by week*/
    /*for(int i=0;i<52;i++)
        save_fscanf(f_mon_pop,"%d %d %d %d %d",&mon_pop[i*5],&mon_pop[i*5+1],&mon_pop[i*5+2],&mon_pop[i*5+3],&mon_pop[i*5+4]);*/


    /*opens the file with the number of positive samples for that strain and season*/
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

    /*********************************************************************************************************************************************************************************
    END of Initialisation bit  LOADING DATA
    **********************************************************************************************************************************************************************************/

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

    one_year_SEIR_with_vaccination(result, pop_vec, curr_init_inf, current_state.time_latent, current_state.time_infectious, current_state.parameters.susceptibility, current_contact_regular, current_state.parameters.transmissibility, vaccine_calendar );

    /*transforms the 21 classes dailys epidemics in weekly 5 AG ones to match RCGP data*/
    days_to_weeks_5AG(result,result_by_week);

    /*curr_psi=0.00001;*/
    lv = log_likelihood_hyper_poisson(current_state.parameters.epsilon, current_state.parameters.psi, result_by_week, ili, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

    auto proposal_state = proposal::initialize( 9 );

    /********************************************************************************************************************************************************
    END initialisation point to start the MCMC
    *********************************************************************************************************************************************************/

    /**************************************************************************************************************************************************************
    Start of the MCMC
    **************************************************************************************************************************************************************/

    step_mat=1;            /*number of contacts exchanged*/
    double p_ac_mat=0.10;          /*prob to redraw matrices*/

    /*mcmc_chain_length=10000;
    acceptance=234;
    freq_sampling=100;
    burn_in=-1;*/

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

        /*proposal_haario(current_state.parameters,proposed_par,chol_emp_cov,chol_ini,100,0.05);*/
        auto proposed_par = proposal::haario_adapt_scale(
                current_state.parameters,
                proposal_state.chol_emp_cov,
                proposal_state.chol_ini,100,0.05, 
                proposal_state.adaptive_scaling, r );

        /*translate into an initial infected population*/
        for(int i=0;i<NAG;i++)
            prop_init_inf[i]=pow(10,proposed_par.init_pop);

        auto prop_c = curr_c;

        /*do swap of contacts step_mat times (reduce or increase to change 'distance' of new matrix from current)*/
        if(R::runif(0,1) < p_ac_mat)
            prop_c = contacts::bootstrap_contacts( std::move(prop_c),
                    polymod_data, step_mat, r );

        auto prop_contact_regular = 
            contacts::to_symmetric_matrix( prop_c, age_data );

        /*one_year_SEIR_without_vaccination(result, pop_vec, prop_init_inf, prop_current_state.time_latent, prop_ti, prop_s_profile, prop_contact_regular, prop_q);*/
        one_year_SEIR_with_vaccination(result, pop_vec, prop_init_inf, current_state.time_latent, current_state.time_infectious, proposed_par.susceptibility, prop_contact_regular, proposed_par.transmissibility, vaccine_calendar );

        /*transforms the 21 classes dailys epidemics in weekly 5 AG ones to match RCGP data*/
        days_to_weeks_5AG(result,result_by_week);

        /*computes the associated likelihood with the proposed values*/
        prop_likelihood=log_likelihood_hyper_poisson(proposed_par.epsilon, proposed_par.psi, result_by_week, ili, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

        /*Acceptance rate include the likelihood and the prior but no correction for the proposal as we use a symmetrical RW*/
        // Make sure accept works with -inf prior
        // MCMC-R alternative prior?
        /*my_acceptance_rate=exp(prop_likelihood-lv+
                log_prior(proposed_par, current_state.parameters, vm.count("prior-susceptibility") ));*/
        my_acceptance_rate=exp(prop_likelihood-lv+
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
            lv=prop_likelihood;

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
    return results;
}
