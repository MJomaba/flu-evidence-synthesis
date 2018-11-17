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

            Eigen::MatrixXd cholesky_I; // Cholesky decomposition of identity matrix

            double adaptive_scaling;
            double conv_scaling;

            size_t no_accepted;
            size_t no_adaptive;
            bool adaptive_step;

            double m;
            double delta;
            double lambda;

            size_t initial_weight = 0; // Weight of initial covariance matrix, used when a correlation matrix is passed along

            proposal_state_t() {
                adaptive_scaling = 0.3;
                conv_scaling = 0.001;

                no_accepted = 0;
                no_adaptive = 0;
                adaptive_step = false;
            }
        };

        proposal_state_t initialize( size_t dim );

        proposal_state_t initialize(Eigen::VectorXd means, Eigen::MatrixXd cov, size_t weight);
        
        proposal_state_t update( proposal_state_t&& state,
                const Eigen::VectorXd &parameters,
                int k );

        /// Was the latest mcmc sample accepted or not
        proposal_state_t accepted( proposal_state_t&& state, 
                bool accepted, int k );

        proposal_state_t update( proposal_state_t&& state,
                const parameter_set &parameters,
                int k );

        /// The original haario version
        Eigen::VectorXd haario( size_t k, 
                const Eigen::VectorXd &current, 
                const Eigen::MatrixXd &chol_de, 
                double epsilon );

        Eigen::VectorXd haario_adapt_scale( const Eigen::VectorXd &current, 
                const Eigen::MatrixXd &chol_de, 
                const Eigen::MatrixXd &chol_ini, 
                double beta, double adapt_scale );

        /// The sherlock 2010 algorithm 6B
        Eigen::VectorXd sherlock( size_t k, 
                const Eigen::VectorXd &current, 
                proposal_state_t &state );
    }
}
#endif

 
