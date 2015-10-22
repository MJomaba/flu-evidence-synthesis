#ifndef FLU_PROPOSAL_HH
#define FLU_PROPOSAL_HH

#include "rcppwrap.h"
#include<RcppEigen.h>

#include <boost/numeric/ublas/matrix.hpp>

#include "model.h"
namespace flu {
    /**
     * \brief Functions to keep track of proposal distribution
     */
    namespace proposal {
        Eigen::VectorXd updateMeans( const Eigen::VectorXd &means,
                const Eigen::VectorXd &v, size_t n );

        /**
         * Incrementally update the covariance matrix
         *
         * This follows the method in:
         * https://root.cern.ch/root/html/TPrincipal.html under AddRow
         * i.e. refer to CERN 72-21 pp. 54-106
         */
        Eigen::MatrixXd updateCovariance( const Eigen::MatrixXd &cov, 
                const Eigen::VectorXd &v, 
                const Eigen::VectorXd &means, 
                size_t n );

        namespace bu = boost::numeric::ublas;
        struct proposal_state_t
        {
            //! Keep track of means of the mcmc samples
            Eigen::VectorXd means_parameters;
            Eigen::MatrixXd emp_cov_matrix, chol_emp_cov, chol_ini;

            double adaptive_scaling = 0.3;
            double past_acceptance = 234;
            double conv_scaling = 0.001;
            int acceptance = 234;
        };

        proposal_state_t initialize( size_t dim );
        
        proposal_state_t update( proposal_state_t&& state,
                const Eigen::VectorXd &parameters,
                int k );

        proposal_state_t update( proposal_state_t&& state,
                const parameter_set &parameters,
                int k );

        Eigen::VectorXd haario_adapt_scale( const Eigen::VectorXd &current, 
                const Eigen::MatrixXd &chol_de, 
                const Eigen::MatrixXd &chol_ini, 
                int n, double beta, double adapt_scale );

        parameter_set haario_adapt_scale( const parameter_set &current, 
                const Eigen::MatrixXd &chol_de, 
                const Eigen::MatrixXd &chol_ini, 
                int n, double beta, double adapt_scale );

    }
}
#endif

 
