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

mcmc_result_inference_t inference_cppWithProposal( std::vector<size_t> demography,
        std::vector<size_t> age_group_limits,
        Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
        Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
        flu::vaccine::vaccine_t vaccine_calendar,
        Eigen::MatrixXi polymod_data,
        const Eigen::VectorXd &initial, const Eigen::VectorXd &initial_contact_ids, 
        proposal::proposal_state_t proposal_state,
        Eigen::MatrixXd mapping,
        Eigen::VectorXd risk_ratios,
        Eigen::VectorXd epsilon_index,
        size_t psi_index,
        size_t transmissibility_index,
        Eigen::VectorXd susceptibility_index,
        size_t initial_infected_index,
        Rcpp::Function lprior,
        bool pass_prior,
        Rcpp::Function lpeak_prior,
        bool pass_peak,
        size_t no_age_groups,
        size_t no_risk_groups,
        bool uk_prior,
        size_t nburn = 0,
        size_t nbatch = 1000, size_t blen = 1, size_t depth = 3)
{
    //Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> mapping,
    mcmc_result_inference_t results;
    results.batch = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>( nbatch, initial.size() );
    results.llikelihoods = Eigen::VectorXd( nbatch );
    results.contact_ids = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>( nbatch, polymod_data.rows() );

    int step_mat;
    double prop_likelihood;

    double my_acceptance_rate;

    flu::data::age_data_t age_data;
    age_data.age_sizes = demography;
    age_data.age_group_sizes = flu::data::group_age_data( demography,
            age_group_limits );

    auto pop_vec = flu::data::stratify_by_risk( 
            age_data.age_group_sizes, risk_ratios, no_risk_groups);

    /*pop RCGP*/
    //std::vector<double> pop_RCGP = std::vector<double>(5);
    Eigen::VectorXd pop_RCGP = Eigen::VectorXd::Zero(ili.cols());
    for (int i = 0; i < mapping.rows(); ++i) {
      // to_j = weight*from_i
      pop_RCGP((size_t) mapping(i,1)) += mapping(i,2) * pop_vec((size_t) mapping(i,0)); 
    }


    /********************************************************************************************************************************************************
    initialisation point to start the MCMC
    *********************************************************************************************************************************************************/

    auto pars_to_epsilon = [&epsilon_index]( const Eigen::VectorXd &pars )
    {
        Eigen::VectorXd eps(epsilon_index.size());
        for (size_t i = 0; i < epsilon_index.size(); ++i)
            eps[i] = pars[epsilon_index[i]];
        return eps;
    };

    auto pars_to_susceptibility = [&susceptibility_index]( const Eigen::VectorXd &pars )
    {
        Eigen::VectorXd susc(susceptibility_index.size());
        for (size_t i = 0; i < susceptibility_index.size(); ++i)
            susc[i] = pars[susceptibility_index[i]];
        return susc;
    };
    
    /*translate into an initial infected population*/
    std::vector<size_t> contact_ids;
    for (size_t i = 0; i < (size_t)polymod_data.rows(); ++i)
        contact_ids.push_back(initial_contact_ids[i]);

    auto curr_parameters = initial;
    Eigen::VectorXd curr_init_inf = flu::data::stratify_by_risk(
            Eigen::VectorXd::Constant(no_age_groups, pow(10,curr_parameters[initial_infected_index]) ),
            risk_ratios, no_risk_groups);

    if (no_risk_groups < 3)
    {
        pop_vec.conservativeResize(no_age_groups*3);
        curr_init_inf.conservativeResize(no_age_groups*3);
        for( size_t i = no_age_groups*no_risk_groups; i<pop_vec.size(); ++i)
        {
            pop_vec[i] = 0;
            curr_init_inf[i] = 0;
        }
    } else if (no_risk_groups > 3)
        ::Rf_error("Maximum of three risk groups supported");

    auto polymod = flu::contacts::table_to_contacts( polymod_data, 
            age_group_limits ); 

    auto curr_c = contacts::shuffle_by_id( polymod, 
            contact_ids );

    auto current_contact_regular = 
        contacts::to_symmetric_matrix( curr_c, age_data );

    auto time_latent = 0.8;
    auto time_infectious = 1.8;

    auto result = infectionODE(pop_vec, 
            curr_init_inf,
            time_latent, time_infectious, 
            pars_to_susceptibility(curr_parameters),
            current_contact_regular, curr_parameters[transmissibility_index], 
            vaccine_calendar, 7*24 );
    /*curr_psi=0.00001;*/
    auto d_app = depth;
    auto curr_llikelihood = log_likelihood_hyper_poisson(
            pars_to_epsilon(curr_parameters),
            curr_parameters[psi_index], 
            days_to_weeks_AG(result, mapping, pop_RCGP.size()), 
            ili, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

    double curr_prior = 0;
    double prop_prior = 0;
    auto Rlprior = [&lprior]( const Eigen::VectorXd &pars ) {
        PutRNGstate();
        double lPrior = Rcpp::as<double>(lprior( pars ));
        GetRNGstate();
        return lPrior;
    };
    
    auto Rlpeak_prior = [&lpeak_prior](const boost::posix_time::ptime &time,
                                       const double &value) {
        boost::gregorian::date d = time.date();
        Rcpp::Date t = Rcpp::Date(d.month(), d.day(), d.year());
        PutRNGstate();
        double lPrior = Rcpp::as<double>(lpeak_prior(t, value));
        GetRNGstate();
        return lPrior;
    };

    if (pass_peak) {
        size_t id;
        auto value = result.cases.rowwise().sum().maxCoeff(&id);
        curr_llikelihood += Rlpeak_prior(result.times[id], value);
    }

    auto log_prior_ratio_f = [pass_prior, &Rlprior, &prop_prior, &curr_prior, uk_prior, &epsilon_index, psi_index, transmissibility_index, &susceptibility_index, 
         initial_infected_index](const Eigen::VectorXd &proposed, const Eigen::VectorXd &current, bool susceptibility) {
             if (uk_prior)
                return log_prior(proposed, current, susceptibility);
             else if (pass_prior) {
                 prop_prior = Rlprior(proposed);
                 return prop_prior - curr_prior;
             } else {
                 prop_prior = 0;
                 // Use flat priors
                 for (auto i = 0; i < epsilon_index.size(); ++i) {
                     auto index = epsilon_index[i];
                     if (proposed[index] < 0 || proposed[index] > 1) {
                         prop_prior = log(0);
                         break;
                     }
                 }
                 if (proposed[psi_index] < 0 || proposed[transmissibility_index] < 0)
                     prop_prior = log(0);
                 for (auto i = 0; i < susceptibility_index.size(); ++i) {
                     auto index = susceptibility_index[i];
                     if (!std::isfinite(prop_prior) || proposed[index] < 0 || proposed[index] > 1) {
                         prop_prior = log(0);
                         break;
                     }
                 }
                 return prop_prior - curr_prior;
             }
         };

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
     
        /*
        if (k>=nburn)
        {
          Rcpp::Rcout << "Adaptive scaling: " << proposal_state.adaptive_scaling << std::endl;
          Rcpp::Rcout << "past_acceptance: " << proposal_state.past_acceptance << std::endl;
          Rcpp::Rcout << "conv_scaling: " << proposal_state.conv_scaling << std::endl;
          Rcpp::Rcout << "Acceptance: " << proposal_state.acceptance << std::endl;
          Rcpp::Rcout << proposal_state.emp_cov_matrix << std::endl << std::endl;
        }
        */
        /*
        auto prop_parameters = proposal::haario_adapt_scale(
                curr_parameters,
                proposal_state.chol_emp_cov,
                proposal_state.chol_ini,0.05, 
                proposal_state.adaptive_scaling );*/

        auto prop_parameters = proposal::sherlock( k,
                curr_parameters,
                proposal_state );

        auto prior_ratio = 
            log_prior_ratio_f(prop_parameters, curr_parameters, false );

        if (!std::isfinite(prior_ratio))
        {
            //Rcpp::Rcout << "Invalid proposed par" << std::endl;
            // TODO: code duplication with failure of acceptance
            proposal_state = proposal::accepted( 
                    std::move(proposal_state), 
                    false, k );
        } else {
            /*translate into an initial infected population*/
            
            Eigen::VectorXd prop_init_inf = flu::data::stratify_by_risk(
                    Eigen::VectorXd::Constant(no_age_groups, pow(10,prop_parameters[initial_infected_index]) ),
                    risk_ratios, no_risk_groups);
            if (no_risk_groups < 3)
            {
                prop_init_inf.conservativeResize(no_age_groups*3);
                for( size_t i = no_age_groups*no_risk_groups; i<prop_init_inf.size(); ++i)
                {
                    prop_init_inf[i] = 0;
                }
            }
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

            result = infectionODE(pop_vec, 
                    prop_init_inf, 
                    time_latent, time_infectious, 
                    pars_to_susceptibility(prop_parameters),
                    prop_contact_regular, prop_parameters[transmissibility_index], 
                    vaccine_calendar, 7*24 );
            
            prop_likelihood = 0;
            if (pass_peak) {
              size_t id;
              auto value = result.cases.rowwise().sum().maxCoeff(&id);
              prop_likelihood = Rlpeak_prior(result.times[id], value);
            }

            /*computes the associated likelihood with the proposed values*/
            prop_likelihood += log_likelihood_hyper_poisson(
                    pars_to_epsilon(prop_parameters), 
                    prop_parameters[psi_index], 
                    days_to_weeks_AG(result, mapping, pop_RCGP.size()), 
                    ili, mon_pop, n_pos, n_samples, pop_RCGP, d_app);

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
                proposal_state = proposal::accepted( 
                        std::move(proposal_state), true, k );

                curr_prior = prop_prior;
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
                proposal_state = proposal::accepted( 
                        std::move(proposal_state), false, k );
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

// [[Rcpp::export(name=".inference_cpp_with_covariance")]]
mcmc_result_inference_t inference_cppWithCovariance( std::vector<size_t> demography, std::vector<size_t> age_group_limits,
        Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
        Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
        flu::vaccine::vaccine_t vaccine_calendar,
        Eigen::MatrixXi polymod_data,
        Eigen::VectorXd initial_parameters, Eigen::VectorXd initial_contact_ids,
        Eigen::VectorXd means, Eigen::MatrixXd covariance, size_t covariance_weight,
        Eigen::MatrixXd mapping,
        Eigen::VectorXd risk_ratios,
        Eigen::VectorXd epsilon_index,
        size_t psi_index,
        size_t transmissibility_index,
        Eigen::VectorXd susceptibility_index,
        size_t initial_infected_index,
        Rcpp::Function lprior,
        bool pass_prior,
        Rcpp::Function lpeak_prior,
        bool pass_peak,
        size_t no_age_groups,
        size_t no_risk_groups,
        bool uk_prior,
        size_t nburn = 0,
        size_t nbatch = 1000, size_t blen = 1, size_t depth = 3)
{
    auto proposal_state = proposal::initialize(means, covariance, covariance_weight);
    return inference_cppWithProposal( demography,
        age_group_limits,
        ili, mon_pop, 
        n_pos, n_samples, 
        vaccine_calendar,
        polymod_data,
        initial_parameters, initial_contact_ids, proposal_state,
        mapping,
        risk_ratios,
        epsilon_index,
        psi_index,
        transmissibility_index,
        susceptibility_index,
        initial_infected_index,
        lprior,
        pass_prior,
        lpeak_prior,
        pass_peak,
        no_age_groups,
        no_risk_groups,
        uk_prior,
        nburn,
        nbatch, blen, depth);
}

//' MCMC based inference of the parameter values given the different data sets
//'
//' @param demography A vector with the population size by each age {0,1,..}
//' @param ili The number of Influenza-like illness cases per week
//' @param mon_pop The number of people monitored for ili
//' @param n_pos The number of positive samples for the given strain (per week)
//' @param n_samples The total number of samples tested 
//' @param vaccine_calendar A vaccine calendar valid for that year
//' @param polymod_data Contact data for different age groups
//' @param initial Vector with starting parameter values
//' @param mapping Group mapping from model groups to data groups
//' @param risk_ratios Risk ratios to convert to and from population groups
//' @param no_age_groups Number of age groups
//' @param no_risk_groups Number of risk groups
//' @param mapping Group mapping from model groups to data groups
//' @param nburn Number of iterations of burn in
//' @param nbatch Number of batches to run (number of samples to return)
//' @param blen Length of each batch
//' 
//' @return Returns a list with the accepted samples and the corresponding llikelihood values and a matrix (contact.ids) containing the ids (row number) of the contacts data used to build the contact matrix.
// [[Rcpp::export(name=".inference_cpp")]]
mcmc_result_inference_t inference_cpp( std::vector<size_t> demography, std::vector<size_t> age_group_limits,
        Eigen::MatrixXi ili, Eigen::MatrixXi mon_pop, 
        Eigen::MatrixXi n_pos, Eigen::MatrixXi n_samples, 
        flu::vaccine::vaccine_t vaccine_calendar,
        Eigen::MatrixXi polymod_data,
        Eigen::VectorXd initial,
        Eigen::MatrixXd mapping,
        Eigen::VectorXd risk_ratios,
        Eigen::VectorXd epsilon_index,
        size_t psi_index,
        size_t transmissibility_index,
        Eigen::VectorXd susceptibility_index,
        size_t initial_infected_index,
        Rcpp::Function lprior,
        bool pass_prior,
        Rcpp::Function lpeak_prior,
        bool pass_peak,
        size_t no_age_groups,
        size_t no_risk_groups,
        bool uk_prior,
        size_t nburn = 0,
        size_t nbatch = 1000, size_t blen = 1, size_t depth = 3)
{
    auto proposal_state = proposal::initialize( initial.size() );
    Eigen::VectorXd contact_ids(polymod_data.rows());
    for (size_t i = 0; i < (size_t)polymod_data.rows(); ++i)
        contact_ids[i] = i+1;
 
  
    return inference_cppWithProposal( demography,
            age_group_limits,
            ili, mon_pop, 
            n_pos, n_samples, 
            vaccine_calendar,
            polymod_data,
            initial, contact_ids, proposal_state,
            mapping,
            risk_ratios,
            epsilon_index,
            psi_index,
            transmissibility_index,
            susceptibility_index,
            initial_infected_index,
            lprior,
            pass_prior,
            lpeak_prior,
            pass_peak,
            no_age_groups,
            no_risk_groups,
            uk_prior,
            nburn,
            nbatch, blen, depth);
}


double dmultinomial( const Eigen::VectorXi &x, int size, 
        const Eigen::VectorXd &prob, 
        bool use_log = false )
{
    double loglik = 0.0;

    loglik = R::lgammafn(size+1);

    for (size_t i=0; i < x.size(); ++i)
        loglik += x[i]*log(prob[i]) - R::lgammafn(x[i]+1);
    if (use_log)
        return loglik;
    return exp(loglik);
}

//' Probability density function for multinomial distribution
//'
//' @param x The counts
//' @param size The total size from which is being samples
//' @param prob Probabilities of each different outcome
//' @param use_log Whether to return logarithm probability
//'
//' @return The probability of getting the counts, given the total size and probability of drawing each.
//'
// [[Rcpp::export(name="dmultinom.cpp")]]
double dmultinomialCPP( Eigen::VectorXi x, int size, Eigen::VectorXd prob, 
        bool use_log = false )
{
    return dmultinomial( x, size, prob, use_log );
}
