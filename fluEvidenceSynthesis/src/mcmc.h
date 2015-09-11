#ifndef FLU_MCMC_HH
#define FLU_MCMC_HH

#include "proposal.h"

struct mcmc_result_t
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        batch;
    Eigen::VectorXd llikelihoods;
};

template<typename Func>
mcmc_result_t adaptiveMCMC( const Func &lprior, const Func &llikelihood,
        size_t nburn,
        const Eigen::VectorXd &initial, 
        size_t nbatch, size_t blen = 1 )
{
    mcmc_result_t result;
    result.batch = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>( nbatch, initial.size() );
    result.llikelihoods = Eigen::VectorXd( nbatch );

    auto curr_parameters = initial;
    auto curr_lprior = lprior( curr_parameters );
    auto curr_llikelihood = llikelihood( curr_parameters );
    auto adapt_rate = 100;
    auto proposal_state = proposal::initialize( initial.size() );

    auto mcmc_chain_length = nbatch*blen;
    for(int k=1; k<=mcmc_chain_length + nburn; k++)
    {

        /*update of the variance-covariance matrix and the mean vector*/
        proposal_state = proposal::update( std::move( proposal_state ),
                curr_parameters, k );

        auto prop_parameters = proposal::haario_adapt_scale(
                curr_parameters,
                proposal_state.chol_emp_cov,
                proposal_state.chol_ini,100,0.05, 
                proposal_state.adaptive_scaling );

        auto prop_lprior = 
            lprior(curr_parameters);

        auto prop_llikelihood = llikelihood( prop_parameters );

        auto my_acceptance_rate = 0.0;
        if (std::isinf(prop_llikelihood) && std::isinf(curr_llikelihood) )
            my_acceptance_rate = exp(prop_lprior-curr_lprior); // We want to explore and find a non infinite likelihood
        else 
            my_acceptance_rate=
                exp(prop_llikelihood-curr_llikelihood+
                        prop_lprior - curr_lprior);

        if(R::runif(0,1)<my_acceptance_rate) /*with prior*/
        {
            /*update the acceptance rate*/
            proposal_state.acceptance++;
            if(k>=adapt_rate)
                proposal_state.adaptive_scaling
                    += 0.766*proposal_state.conv_scaling;

            curr_parameters = prop_parameters;

            /*update current likelihood*/
            curr_llikelihood=prop_llikelihood;
            curr_lprior=prop_lprior;
        }
        else /*if reject*/
        {
            if(k>=adapt_rate)
                proposal_state.adaptive_scaling
                    -=0.234*proposal_state.conv_scaling;
        }
        if(k%blen==0 && k>=nburn)
        {
            // Add results
            size_t i = (k-nburn)/blen;
            result.llikelihoods[i] = curr_llikelihood;
            result.batch.row( i ) = curr_parameters;
        }
    }
    return result;
}

#endif
