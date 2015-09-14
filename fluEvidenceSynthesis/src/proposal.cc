#include "proposal.h"
#include <boost/numeric/ublas/matrix.hpp>
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
            // Specialized case
            /*if (init_cov_matrix.rows()==9)
            {
                init_cov_matrix(3,3) = 1e-16;
                init_cov_matrix(4,4) = 0.000001;
                init_cov_matrix(5,5) = 0.000007;
                init_cov_matrix(6,6) = 0.000007;
                init_cov_matrix(7,7) = 0.000007;
                init_cov_matrix(8,8) = 0.00003;
            }*/

            state.chol_ini = Eigen::LLT<Eigen::MatrixXd>(
                    init_cov_matrix).matrixL();
            state.chol_emp_cov = state.chol_ini;

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
            if(k%1000==0)
            {
                state.chol_emp_cov = Eigen::LLT<Eigen::MatrixXd>(
                        state.emp_cov_matrix).matrixL();

                state.past_acceptance=state.acceptance;
                state.acceptance=0;

                state.conv_scaling/=1.005;
            }
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

         Eigen::VectorXd haario_adapt_scale( const Eigen::VectorXd &current, 
                const Eigen::MatrixXd &chol_de, 
                const Eigen::MatrixXd &chol_ini, 
                int n, double beta, double adapt_scale )
        {
            auto proposed = Eigen::VectorXd( current.size() );
            auto normal_draw = Eigen::VectorXd( current.size() );
            auto normal_add_draw = Eigen::VectorXd( current.size() );
            auto correlated_draw = Eigen::VectorXd( current.size() );
            auto correlated_fix = Eigen::VectorXd( current.size() );
            double unif1, unif2;
            double un_moins_beta;
            un_moins_beta=1-beta;

            // TODO create random multivariate draw and use with
            // both chol_de and chol_ini
            /*drawing of the needed N(0,1) samples using Box-Muller*/
            for(size_t i=0;i<current.size();i++)
            {
                unif1=R::runif(0,1);
                unif2=R::runif(0,1);
                normal_draw[i]=adapt_scale*2.38/3*sqrt(-2*log(unif1))*sin(twopi*unif2); /*3 = sqrt(9)*/
                /*drawing of the needed N(0,1) samples using Box-Muller*/
                normal_add_draw[i]=sqrt(-2*log(unif1))*cos(twopi*unif2);
            }

            /*transforming the numbers generated with the Cholesky mat to get the correlated samples*/
            // TODO use vector::zero
            for(size_t i=0;i<current.size();i++)
                correlated_draw[i]=0;

            for(size_t i=0;i<current.size();i++)
                for(size_t j=0;j<=i;j++)
                    correlated_draw[i]+=chol_de(i,j)*normal_draw[j];

            // TODO use vector::zero
            for(size_t i=0;i<current.size();i++)
                correlated_fix[i]=0;

            for(size_t i=0;i<current.size();i++)
                for(size_t j=0;j<=i;j++)
                    correlated_fix[i]+=chol_ini(i,j)*normal_add_draw[j];

            /*new proposed values*/
            proposed = current + un_moins_beta*correlated_draw + beta*correlated_fix;

            return proposed;
        }


        parameter_set haario_adapt_scale( const parameter_set &current, 
                const Eigen::MatrixXd &chol_de, 
                const Eigen::MatrixXd &chol_ini, 
                int n, double beta, double adapt_scale )
        {
            parameter_set proposed;
            double normal_draw[9];
            double normal_add_draw[9], correlated_draw[9], correlated_fix[9];
            double unif1, unif2;
            int i, j;
            double un_moins_beta;

            un_moins_beta=1-beta;

            // TODO create random multivariate draw and use with
            // both chol_de and chol_ini
            /*drawing of the 9 N(0,1) samples using Box-Muller*/
            for(i=0;i<4;i++)
            {
                unif1=R::runif(0,1);
                unif2=R::runif(0,1);
                normal_draw[i*2]=adapt_scale*2.38/3*sqrt(-2*log(unif1))*sin(twopi*unif2); /*3 = sqrt(9)*/
                normal_draw[i*2+1]=adapt_scale*2.38/3*sqrt(-2*log(unif1))*cos(twopi*unif2);
            }

            unif1=R::runif(0,1);
            unif2=R::runif(0,1);
            normal_draw[8]=adapt_scale*2.38/3*sqrt(-2*log(unif1))*sin(twopi*unif2);
            normal_add_draw[8]=sqrt(-2*log(unif1))*cos(twopi*unif2);

            /*drawing of the 9 N(0,1) samples using Box-Muller*/
            for(i=0;i<4;i++)
            {
                unif1=R::runif(0,1);
                unif2=R::runif(0,1);
                normal_add_draw[i*2]=sqrt(-2*log(unif1))*sin(twopi*unif2);
                normal_add_draw[i*2+1]=sqrt(-2*log(unif1))*cos(twopi*unif2);
            }

            /*transforming the numbers generated with the Cholesky mat to get the correlated samples*/
            for(i=0;i<9;i++)
                correlated_draw[i]=0;

            for(i=0;i<9;i++)
                for(j=0;j<=i;j++)
                    correlated_draw[i]+=chol_de(i,j)*normal_draw[j];

            for(i=0;i<9;i++)
                correlated_fix[i]=0;

            for(i=0;i<9;i++)
                for(j=0;j<=i;j++)
                    correlated_fix[i]+=chol_ini(i,j)*normal_add_draw[j];

            /*new proposed values*/
            proposed.epsilon[0]=current.epsilon[0]+un_moins_beta*correlated_draw[0]+beta*correlated_fix[0];
            proposed.epsilon[1]=proposed.epsilon[0];
            proposed.epsilon[2]=current.epsilon[2]+un_moins_beta*correlated_draw[1]+beta*correlated_fix[1];
            proposed.epsilon[3]=proposed.epsilon[2];
            proposed.epsilon[4]=current.epsilon[4]+un_moins_beta*correlated_draw[2]+beta*correlated_fix[2];

            proposed.psi=current.psi+un_moins_beta*correlated_draw[3]+beta*correlated_fix[3];

            proposed.transmissibility=current.transmissibility+un_moins_beta*correlated_draw[4]+beta*correlated_fix[4];

            proposed.susceptibility[0]=current.susceptibility[0]+un_moins_beta*correlated_draw[5]+beta*correlated_fix[5];
            proposed.susceptibility[1]=proposed.susceptibility[0];
            proposed.susceptibility[2]=proposed.susceptibility[0];
            proposed.susceptibility[3]=current.susceptibility[3]+un_moins_beta*correlated_draw[6]+beta*correlated_fix[6];
            proposed.susceptibility[4]=proposed.susceptibility[3];
            proposed.susceptibility[5]=proposed.susceptibility[3];
            proposed.susceptibility[6]=current.susceptibility[6]+un_moins_beta*correlated_draw[7]+beta*correlated_fix[7];

            proposed.init_pop=current.init_pop+un_moins_beta*correlated_draw[8]+beta*normal_add_draw[8];

            return proposed;
        }
    }
}

