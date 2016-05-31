#include "proposal.h"
#include <boost/numeric/ublas/matrix.hpp>

#define twopi 6.283185

namespace flu {
    /**
     * \brief Functions to keep track of proposal distribution
     */
    namespace proposal {
        Eigen::VectorXd updateMeans( const Eigen::VectorXd &means,
                const Eigen::VectorXd &v, size_t n )
        {
            if (n==1)
                return v;
            return means + 1.0/n*(v-means);
        }

        Eigen::MatrixXd updateCovariance( const Eigen::MatrixXd &cov, 
                const Eigen::VectorXd &v, 
                const Eigen::VectorXd &means, 
                size_t n )
        {
            if (n==1)
                return Eigen::MatrixXd::Zero( v.size(), v.size() );
            return cov + (1.0/(n-1.0))*((v-means)*((v-means).transpose())) 
                - (1.0/n)*cov;
        }
        
        proposal_state_t initialize( size_t dim )
        {
            proposal_state_t state;
            state.means_parameters = Eigen::VectorXd::Zero( dim );

            Eigen::MatrixXd init_cov_matrix = 
                Eigen::MatrixXd::Zero(dim,dim);

            state.emp_cov_matrix =
                Eigen::MatrixXd::Zero(dim,dim);

            state.chol_emp_cov.resize( dim, dim );
            state.chol_ini.resize( dim, dim );

            for(int i=0; i<init_cov_matrix.rows(); i++)
            {
                init_cov_matrix(i,i) = 0.0000001;
            }

            state.chol_ini = Eigen::LLT<Eigen::MatrixXd>(
                    init_cov_matrix).matrixL();
            state.chol_emp_cov = state.chol_ini;

            state.cholesky_I = Eigen::LLT<Eigen::MatrixXd>(
                    Eigen::MatrixXd::Identity(dim,dim))
                .matrixL();

            state.m = 2.38/sqrt(dim);
            state.delta = state.m/100;
            state.lambda = 0.001/sqrt(dim);

            return state;
        }

        /// Was the latest mcmc sample accepted or not
        proposal_state_t accepted( proposal_state_t&& state, 
                bool accepted, int k )
        {
            if (state.adaptive_step)
                ++state.no_adaptive;

            if (accepted)
            {
                ++state.no_accepted;
                if (k>=100)
                    state.adaptive_scaling
                        += 0.766*state.conv_scaling;
                if (state.adaptive_step)
                    state.m += 2.3*state.delta/sqrt(state.no_adaptive);
            }
            else {
                if (state.adaptive_step)
                    state.m -= state.delta/sqrt(state.no_adaptive);
            
                if (k>=100)
                    state.adaptive_scaling
                        -= 0.234*state.conv_scaling;
            } 

            return state;
        }

        proposal_state_t update( proposal_state_t&& state,
                const Eigen::VectorXd &parameters,
                int k )
        {
            /*update of the variance-covariance matrix and the mean vector*/
            state.means_parameters = updateMeans( 
                    state.means_parameters, parameters, k );
            state.emp_cov_matrix = updateCovariance( state.emp_cov_matrix,
                    parameters, state.means_parameters, k );
            /*adjust variance for MCMC parameters*/
            /*if(k%100==0)
            {*/
                state.chol_emp_cov = Eigen::LLT<Eigen::MatrixXd>(
                        state.emp_cov_matrix).matrixL();

                state.conv_scaling/=1.005;
            //}
            return state;
        }

        proposal_state_t update( proposal_state_t&& state,
                const parameter_set &parameters,
                int k ) 
        {
            // Vectorize parameters
            //vectorize parameters
            Eigen::VectorXd pars_v( 9 );
            pars_v << parameters.epsilon[0], parameters.epsilon[2],
                   parameters.epsilon[4],
                   parameters.psi, parameters.transmissibility,
                   parameters.susceptibility[0], 
                   parameters.susceptibility[3],
                   parameters.susceptibility[6], parameters.init_pop;

            return update( std::move(state), pars_v, k );
        }

        Eigen::VectorXd haario( size_t k,
                const Eigen::VectorXd &current, 
                const Eigen::MatrixXd &chol_de, 
                double epsilon )
        {
            auto proposed = Eigen::VectorXd( current.size() );
            auto normal_draw = Eigen::VectorXd( current.size() );
            auto normal_add_draw = Eigen::VectorXd( current.size() );
            Eigen::VectorXd correlated_draw = 
                Eigen::VectorXd::Zero( current.size() );
            double unif1, unif2;

            auto sqrtd = sqrt( current.size() );
            auto sd = 2.38/sqrtd;
            auto esd = epsilon*sd;

            // TODO create random multivariate draw and use with
            // both chol_de and chol_ini
            /*drawing of the needed N(0,1) samples using Box-Muller*/
            for(int i=0;i<current.size();i++)
            {
                unif1=R::runif(0,1);
                unif2=R::runif(0,1);
                normal_draw[i]=sqrt(-2*log(unif1))*sin(twopi*unif2); /*3 = sqrt(9)*/
                /*drawing of the needed N(0,1) samples using Box-Muller*/
                normal_add_draw[i]=sqrt(-2*log(unif1))*cos(twopi*unif2);
            }

            if (k<1000)
            {
                proposed = current + esd*normal_add_draw;
            }

            /*transforming the numbers generated with the Cholesky mat to get the correlated samples*/
            for(int i=0;i<current.size();i++)
                for(int j=0;j<=i;j++)
                    correlated_draw[i]+=chol_de(i,j)*normal_draw[j];

            // Identity matrix should mean that Id*normal_add_draw = normal_add_draw
            /*new proposed values*/
            proposed = current + sd*correlated_draw + esd*normal_add_draw;
            return proposed;
        }





         Eigen::VectorXd haario_adapt_scale( const Eigen::VectorXd &current, 
                const Eigen::MatrixXd &chol_de, 
                const Eigen::MatrixXd &chol_ini, 
                double beta, double adapt_scale )
        {
            auto proposed = Eigen::VectorXd( current.size() );
            auto normal_draw = Eigen::VectorXd( current.size() );
            auto normal_add_draw = Eigen::VectorXd( current.size() );
            Eigen::VectorXd correlated_draw = 
                Eigen::VectorXd::Zero( current.size() );
            Eigen::VectorXd correlated_fix = 
                Eigen::VectorXd::Zero( current.size() );
            double unif1, unif2;
            double un_moins_beta;
            un_moins_beta=1-beta;

            auto sqrtd = sqrt( current.size() );

            // TODO create random multivariate draw and use with
            // both chol_de and chol_ini
            /*drawing of the needed N(0,1) samples using Box-Muller*/
            for(int i=0;i<current.size();i++)
            {
                unif1=R::runif(0,1);
                unif2=R::runif(0,1);
                normal_draw[i]=adapt_scale*2.38/sqrtd*sqrt(-2*log(unif1))*sin(twopi*unif2); /*3 = sqrt(9)*/
                /*drawing of the needed N(0,1) samples using Box-Muller*/
                normal_add_draw[i]=sqrt(-2*log(unif1))*cos(twopi*unif2);
            }

            /*transforming the numbers generated with the Cholesky mat to get the correlated samples*/
            for(int i=0;i<current.size();i++)
                for(int j=0;j<=i;j++)
                    correlated_draw[i]+=chol_de(i,j)*normal_draw[j];

            for(int i=0;i<current.size();i++)
                for(int j=0;j<=i;j++)
                    correlated_fix[i]+=chol_ini(i,j)*normal_add_draw[j];

            /*new proposed values*/
            proposed = current + un_moins_beta*correlated_draw + beta*correlated_fix;

            return proposed;
        }

        /// The sherlock 2010 algorithm 6B
        Eigen::VectorXd sherlock( size_t k, 
                const Eigen::VectorXd &current, 
                proposal_state_t &state ) {

            auto proposed = Eigen::VectorXd( current.size() );
            auto normal_draw = Eigen::VectorXd( current.size() );
            Eigen::VectorXd correlated = 
                Eigen::VectorXd::Zero( current.size() );
            //auto normal_add_draw = Eigen::VectorXd( current.size() );
 
            for(int i=0;i<current.size();i++)
            {
                auto unif1=R::runif(0,1);
                auto unif2=R::runif(0,1);
                normal_draw[i]=sqrt(-2*log(unif1))*sin(twopi*unif2);
                /*drawing of the needed N(0,1) samples using Box-Muller*/
                //normal_add_draw[i]=sqrt(-2*log(unif1))*cos(twopi*unif2);
            }

            if (state.no_accepted<100 || R::runif(0,1)<0.05)
            {
                state.adaptive_step = false;
                for(int i=0;i<current.size();i++)
                    for(int j=0;j<=i;j++)
                        correlated[i]+=state.cholesky_I(i,j)*normal_draw[j];

                 proposed = current + state.lambda*correlated;
            } else {
                state.adaptive_step = true;
                for(int i=0;i<current.size();i++)
                    for(int j=0;j<=i;j++)
                        correlated[i]+=state.chol_emp_cov(i,j)*normal_draw[j];

                proposed = current + state.m*correlated;
            }
            return proposed;
        }


    }
}

