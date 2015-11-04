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
        std::vector<double> separate_into_risk_groups( 
                const std::vector<size_t> &group_sizes )
        {
            /*opens the file with the number of positive samples for that strain and season*/
            //auto f_pop_model = read_file( filepath );

            std::vector<double> pop_vec;
            pop_vec.resize( 21, 0 );
            /*load the size of the age groups for the model that year*/
            //double pop_vec[21];
            for( size_t i=0; i < group_sizes.size(); ++i ) 
                pop_vec[i] = group_sizes[i];

            /*high risk*/
            pop_vec[7]=pop_vec[0]*0.021; /*2.1% in the <1 */
            pop_vec[8]=pop_vec[1]*0.055; /*5.5% in the 1-4 */
            pop_vec[9]=pop_vec[2]*0.098; /*9.8% in the 5-14 */
            pop_vec[10]=pop_vec[3]*0.087; /*8.7% in the 15-24 */
            pop_vec[11]=pop_vec[4]*0.092; /*9.2% in the 25-44 */
            pop_vec[12]=pop_vec[5]*0.183; /*18.3% in the 45-64 */
            pop_vec[13]=pop_vec[6]*0.45; /*45% in the 65+ */

            /*pregnant women (pregnant women not considered in this model*/
            pop_vec[14]=0;
            pop_vec[15]=0;
            pop_vec[16]=0;
            pop_vec[17]=0;
            pop_vec[18]=0;
            pop_vec[19]=0;
            pop_vec[20]=0;

            /*low risk*/
            pop_vec[0]-=pop_vec[7];
            pop_vec[1]-=pop_vec[8];
            pop_vec[2]-=pop_vec[9];
            pop_vec[3]-=pop_vec[10]+pop_vec[17];
            pop_vec[4]-=pop_vec[11]+pop_vec[18];
            pop_vec[5]-=pop_vec[12];
            pop_vec[6]-=pop_vec[13];
            //fclose( f_pop_model );

            return pop_vec;
        }

        std::vector<size_t> group_age_data( 
                const std::vector<size_t> &age_sizes,
                const std::vector<size_t> &limits )
        {
            std::vector<size_t> age_group_sizes(limits.size()+1, 0);

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
