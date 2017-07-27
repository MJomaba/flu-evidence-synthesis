#ifndef FLU_DATA_HH
#define FLU_DATA_HH

#include<string>
#include<vector>

#include "state.h"

namespace flu
{
    namespace data
    {
        /** 
         * \brief Separate the groups sizes into risk groups
         *
         * Example groups are: low risk, high risk and pregnant women
         */
        template<class T>
            Eigen::VectorXd separate_into_risk_groups( 
                    const T &age_groups,
                    const Eigen::MatrixXd &risk )
            {
                Eigen::VectorXd pop_vec( age_groups.size() * (1+risk.rows()) );
                // Fill risk groups
                for (int i = 0; i < risk.rows(); ++i)
                {
                    for (int j = 0; j < age_groups.size(); ++j)
                    {
                        pop_vec[ (1+i)*age_groups.size() + j ] =
                            age_groups[j]*risk(i,j);
                    }
                }

                // Fill low risk with left over
                for( int i = 0; i < age_groups.size(); ++i )
                {
                    pop_vec[i] = age_groups[i];
                    for (int j = 0; j < risk.rows(); ++ j)
                    {
                        pop_vec[i] -= pop_vec[(1+j)*age_groups.size() + i];
                    }
                }

                return pop_vec;
            }

        /** 
         * \brief Separate the groups sizes into risk groups
         */
        template<class T>
            Eigen::VectorXd stratify_by_risk( 
                    const T &age_groups, const Eigen::VectorXd &risk, size_t no_risk_groups )
            {
                Eigen::VectorXd pop_vec(age_groups.size() * no_risk_groups);
                for (auto i = 0; i < no_risk_groups; ++i) {
                    for (auto j = 0; j < age_groups.size(); ++j) {
                        auto id = j+age_groups.size()*i;
                        pop_vec[id] = age_groups[j]*risk[id];
                    }
                }
                return pop_vec;
            }
 
        /**
         * \brief Group population size according to age groups
         */
        Eigen::VectorXi group_age_data( 
                const std::vector<size_t> &age_sizes, 
                const std::vector<size_t> &limits );

        struct age_data_t
        {
            /// Population per year/age
            std::vector<size_t> age_sizes;

            /// Population per age group
            Eigen::VectorXi age_group_sizes;
        };
    }
}

#endif
