#include "data.hh"

#include <stdio.h>

#include "io.hh"
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

        std::vector<double> load_contact_regular( const std::string& contacts_filepath, const state_t &state, int *age_sizes, int *AG_sizes )
        {
            int nc, age_part, AG_part;
            double ww[POLY_PART], mij[49], w_norm[7], cij[49], cij_pro;
            int c_age[POLY_PART], c_we[POLY_PART], c_N1[POLY_PART], c_N2[POLY_PART], c_N3[POLY_PART], c_N4[POLY_PART], c_N5[POLY_PART], c_N6[POLY_PART], c_N7[POLY_PART], c_AG[POLY_PART], c_ni[90];
            int curr_age[POLY_PART], curr_we[POLY_PART], curr_N1[POLY_PART], curr_N2[POLY_PART], curr_N3[POLY_PART], curr_N4[POLY_PART], curr_N5[POLY_PART], curr_N6[POLY_PART], curr_N7[POLY_PART], curr_AG[POLY_PART], curr_ni[90], curr_nwe;
            FILE *contacts_PM;
            int c_nwe = 0;
            contacts_PM=read_file( contacts_filepath );
            for(size_t i=0; i<POLY_PART; i++)
            {
                save_fscanf(contacts_PM,"%d %d %d %d %d %d %d %d %d", &c_age[i], &c_we[i], &c_N1[i], &c_N2[i], &c_N3[i], &c_N4[i], &c_N5[i], &c_N6[i], &c_N7[i]);
                age_part=c_age[i];
                c_ni[age_part]++;
                if(c_we[i]>0) c_nwe++;
                c_AG[i]=0;
                if(age_part>0) c_AG[i]++;
                if(age_part>4) c_AG[i]++;
                if(age_part>14) c_AG[i]++;
                if(age_part>24) c_AG[i]++;
                if(age_part>44) c_AG[i]++;
                if(age_part>64) c_AG[i]++;
            }


            for(size_t i=0;i<90;i++)
                curr_ni[i]=0;
            curr_nwe=0;
            for(size_t i=0; i<POLY_PART; i++)
            {
                nc = state.contact_ids[i];

                age_part=c_age[nc];
                curr_ni[age_part]++;
                if(c_we[nc]>0) curr_nwe++;

                curr_age[i]=c_age[nc];
                curr_AG[i]=c_AG[nc];
                curr_we[i]=c_we[nc];
                curr_N1[i]=c_N1[nc];
                curr_N2[i]=c_N2[nc];
                curr_N3[i]=c_N3[nc];
                curr_N4[i]=c_N4[nc];
                curr_N5[i]=c_N5[nc];
                curr_N6[i]=c_N6[nc];
                curr_N7[i]=c_N7[nc];
            }

            /*update of the weights*/
            for(size_t i=0;i<7;i++)
                w_norm[i]=0;

            for(size_t i=0;i<49;i++)
                mij[i]=0;

            for(size_t i=0; i<POLY_PART; i++)
            {
                age_part=curr_age[i];
                AG_part=curr_AG[i];
                if(curr_we[i]==0)
                    ww[i]=(double)age_sizes[age_part]/curr_ni[age_part]*5/(POLY_PART-curr_nwe);
                else
                    ww[i]=(double)age_sizes[age_part]/curr_ni[age_part]*2/curr_nwe;

                w_norm[AG_part]+=ww[i];
                mij[7*AG_part]+=curr_N1[i]*ww[i];
                mij[7*AG_part+1]+=curr_N2[i]*ww[i];
                mij[7*AG_part+2]+=curr_N3[i]*ww[i];
                mij[7*AG_part+3]+=curr_N4[i]*ww[i];
                mij[7*AG_part+4]+=curr_N5[i]*ww[i];
                mij[7*AG_part+5]+=curr_N6[i]*ww[i];
                mij[7*AG_part+6]+=curr_N7[i]*ww[i];
            }

            for(size_t i=0; i<49; i++)
            {
                if(w_norm[i/7]>0)
                    mij[i]/=w_norm[i/7];
                cij[i]=mij[i]/AG_sizes[i%7];
            }

            std::vector<double> contact_regular( NAG2 );
            for(size_t i=0; i<7; i++)
            {
                contact_regular[i*7+i]=cij[i*7+i];
                for(size_t j=0;j<i;j++)
                {
                    cij_pro=(cij[i*7+j]+cij[j*7+i])/2;
                    contact_regular[i*7+j]=cij_pro;
                    contact_regular[j*7+i]=cij_pro;
                }
            }

            fclose( contacts_PM );

            return contact_regular;
        }

        age_data_t load_age_data( const std::string &path )
        {
            age_data_t age_data;
            auto pop_sizes = read_file( path );
            for(size_t i=0; i<7; i++)
                age_data.age_group_sizes[i]=0;
            for(size_t i=0; i<85; i++)
            {
                save_fscanf(pop_sizes,"%d",&age_data.age_sizes[i]);
                if(i==0)
                    age_data.age_group_sizes[0]=age_data.age_sizes[0];
                else
                    if(i<5)
                        age_data.age_group_sizes[1]+=age_data.age_sizes[i];
                    else
                        if(i<15)
                            age_data.age_group_sizes[2]+=age_data.age_sizes[i];
                        else
                            if(i<25)
                                age_data.age_group_sizes[3]+=age_data.age_sizes[i];
                            else
                                if(i<45)
                                    age_data.age_group_sizes[4]+=age_data.age_sizes[i];
                                else
                                    if(i<65)
                                        age_data.age_group_sizes[5]+=age_data.age_sizes[i];
                                    else
                                        age_data.age_group_sizes[6]+=age_data.age_sizes[i];
            }

            /*put the remaining (85+) in the older AG*/
            int aux;
            save_fscanf(pop_sizes,"%d",&aux);
            age_data.age_group_sizes[6]+=aux;

            return age_data;
        }
    };
};
