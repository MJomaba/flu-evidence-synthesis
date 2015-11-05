#include "data.h"

#include "rcppwrap.h"

#include <stdio.h>
#include <fstream>
#include <deque>

#include "state.h"

namespace flu
{
    namespace data
    {

        Eigen::VectorXd separate_into_risk_groups( 
                const Eigen::VectorXi &age_groups,
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

        Eigen::VectorXi group_age_data( 
                const std::vector<size_t> &age_sizes,
                const std::vector<size_t> &limits )
        {
            Eigen::VectorXi age_group_sizes = 
                Eigen::VectorXi::Zero(limits.size()+1);

            size_t current_age = 0;
            size_t group_count = 0;
            std::deque<size_t> group_barriers = 
                std::deque<size_t>( limits.begin(), limits.end() );

            // Iterate over each line
            for( auto &pop_size : age_sizes )
            {
                if (group_barriers.size() != 0 &&
                        current_age >= group_barriers[0] )
                {
                    // Population should be added to the next age group
                    ++group_count;
                    group_barriers.pop_front();
                }
                age_group_sizes[group_count]+=pop_size;
                ++current_age;
            }

            return age_group_sizes;
        }

    }
}
