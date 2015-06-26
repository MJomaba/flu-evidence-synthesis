#include "data.hh"

#include <stdio.h>
#include <fstream>
#include <deque>

#include "state.hh"

namespace flu
{
    namespace data
    {
        std::vector<double> separate_into_risk_groups( 
                const std::array<size_t,7> &group_sizes )
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

    };
};
