#ifndef FLU_PROPOSAL_HH
#define FLU_PROPOSAL_HH

#include <boost/numeric/ublas/matrix.hpp>

#include "model.hh"
namespace flu {
    /**
     * \brief Functions to keep track of proposal distribution
     */
    namespace proposal {
        namespace bu = boost::numeric::ublas;
        struct proposal_state_t
        {
            //! Keep track of means of the mcmc samples
            std::vector<double> sum_mean_param;
            bu::matrix<double> emp_cov_matrix, sum_corr_param_matrix, chol_emp_cov, chol_ini;

            double adaptive_scaling = 0.3;
            double past_acceptance = 234;
            double conv_scaling = 0.001;
            int acceptance = 234;
        };

        bu::matrix<double> cholesky_factorization(
                const bu::matrix<double> &A);
        
        bu::matrix<double> update_sum_corr(bu::matrix<double> &&corr, 
                const parameter_set &par );
        
        proposal_state_t update( proposal_state_t&& state,
                const parameter_set &parameters,
                int k );

        parameter_set haario_adapt_scale(const parameter_set &current, 
                const bu::matrix<double> &chol_de, const bu::matrix<double> &chol_ini, 
                int n, double beta, double adapt_scale, gsl_rng * gen);

    };
};
#endif

 
