#ifndef FLU_MCMC_HH
#define FLU_MCMC_HH

#include<chrono>

#include "proposal.h"

namespace flu {

struct mcmc_result_t
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        batch;
    Eigen::VectorXd llikelihoods;
};

template<typename Func1, typename Func2, typename Func3, typename Func4>
mcmc_result_t adaptiveMCMC( const Func1 &lprior, const Func2 &llikelihood, 
        const Func3 &outfun, const Func4 &acceptfun,
        size_t nburn,
        const Eigen::VectorXd &initial, 
        size_t nbatch, size_t blen = 1, bool verbose = false )
{
    mcmc_result_t result;
    result.batch = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>( nbatch, initial.size() );
    result.llikelihoods = Eigen::VectorXd( nbatch );

    auto curr_parameters = initial;
    PutRNGstate();
    auto curr_lprior = lprior( curr_parameters );
    auto curr_llikelihood = llikelihood( curr_parameters );
    GetRNGstate();
    auto proposal_state = proposal::initialize( initial.size() );


    if (verbose) {
        Rcpp::Rcout << "Initial LPrior\t" << curr_lprior  << std::endl;
        Rcpp::Rcout << "Initial Llikeli\t" << curr_llikelihood << std::endl;
    }


    size_t sampleCount = 0;
    int k = 0;
    while(sampleCount<nbatch)
    {
        ++k;

        //update of the variance-covariance matrix and the mean vector
        proposal_state = proposal::update( std::move( proposal_state ),
                curr_parameters, k );

        /*
        auto prop_parameters = proposal::haario_adapt_scale(
                curr_parameters,
                proposal_state.chol_emp_cov,
                proposal_state.chol_ini,0.05, 
                proposal_state.adaptive_scaling );*/
        /*
        auto epsilon = 0.001;
        if (k>10000)
            epsilon = 0;
        auto prop_parameters = proposal::haario( k,
                curr_parameters,
                proposal_state.chol_emp_cov, epsilon );
                */
        auto prop_parameters = proposal::sherlock( k,
                curr_parameters,
                proposal_state );
        if (verbose) {
            Rcpp::Rcout << "Proposed parameters\t" << prop_parameters.transpose() <<
                std::endl;
            Rcpp::Rcout << "Var\tvalue\ttime (ns)" << std::endl;
        }
        std::chrono::high_resolution_clock::time_point start_time; 

        if (verbose)
            start_time = std::chrono::high_resolution_clock::now();
        PutRNGstate();
        auto prop_lprior = 
            lprior(prop_parameters);
        GetRNGstate();
        if (verbose)
            Rcpp::Rcout << "LPrior\t" << prop_lprior << "\t" <<
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() -
                        start_time).count() << std::endl;
        auto prop_llikelihood = log(0);

        auto my_acceptance_rate = 0.0;
        if (!std::isinf(prop_lprior)) 
        {
            if (verbose)
                start_time = std::chrono::high_resolution_clock::now();
            PutRNGstate();
            prop_llikelihood = llikelihood( prop_parameters );
            GetRNGstate();
            if (verbose)
                Rcpp::Rcout << "Llikeli\t" << prop_llikelihood << "\t" <<
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() -
                            start_time).count() << std::endl;

            if (std::isinf(prop_llikelihood) && std::isinf(curr_llikelihood) )
                my_acceptance_rate = exp(prop_lprior-curr_lprior); // We want to explore and find a non infinite likelihood
            else
                my_acceptance_rate=
                    exp(prop_llikelihood-curr_llikelihood+
                            prop_lprior - curr_lprior);
        } else {
            if (std::isinf(curr_lprior))
                ::Rf_error("Algorithm stuck on infinite prior");
            my_acceptance_rate = -1.0;
        }
        auto rnd = R::runif(0.0, 1.0);
        if (verbose)
            Rcpp::Rcout << "RND: " << rnd << " rate " << my_acceptance_rate << std::endl;
        if(rnd < my_acceptance_rate) //with prior
        {
            //update the acceptance rate
            proposal_state = proposal::accepted( 
                    std::move(proposal_state), true, k );

            curr_parameters = prop_parameters;

            //update current likelihood
            curr_llikelihood=prop_llikelihood;
            curr_lprior=prop_lprior;
            // Call the accept function
            if (verbose)
                start_time = std::chrono::high_resolution_clock::now();
            PutRNGstate();
            acceptfun();
            GetRNGstate();
            if (verbose)
                Rcpp::Rcout << "acceptfun\t" << 0.0/0.0 << "\t" <<
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() -
                            start_time).count() << std::endl;

        }
        else //if reject
        {
            proposal_state = proposal::accepted( 
                    std::move(proposal_state), false, k );
        }

        if(k%blen==0 && k>=nburn)
        {
            // Add results
            result.llikelihoods[sampleCount] = curr_llikelihood;
            result.batch.row( sampleCount ) = curr_parameters;
            if (verbose)
                start_time = std::chrono::high_resolution_clock::now();
            PutRNGstate();
            outfun();
            GetRNGstate();
            if (verbose)
                Rcpp::Rcout << "outfun\t" << 0.0/0.0 << "\t" <<
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() -
                            start_time).count() << std::endl;

            ++sampleCount;
        }
    }
    return result;
}
}
#endif
