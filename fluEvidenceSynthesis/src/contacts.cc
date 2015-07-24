#include "contacts.h"

#include "rcppwrap.h"

#include <stdio.h>
#include <cassert>

#include "state.h"

namespace flu
{
    namespace contacts {
        contacts_t bootstrap_contacts( contacts_t&& bootstrap,
                const contacts_t &original,
                size_t no )
        {
            for(size_t i=0;i<no;i++)
            {
                auto alea1=(size_t) R::runif(0,POLY_PART);
                auto alea2=(size_t) R::runif(0,POLY_PART);

                bootstrap.ni[bootstrap.contacts[alea1].age]--;
                if(bootstrap.contacts[alea1].weekend) bootstrap.nwe--;

                bootstrap.contacts[alea1]=original.contacts[alea2];

                bootstrap.ni[bootstrap.contacts[alea1].age]++;
                if(bootstrap.contacts[alea1].weekend) bootstrap.nwe++;
            }
            return bootstrap;
        }

        contacts_t shuffle_by_id( const contacts_t &sorted_c, const std::vector<size_t> &ids )
        {
            contacts_t shuffled_c;

            /*translate into an initial infected population*/
            for(size_t i=0;i<90;i++)
                shuffled_c.ni[i]=0;
            shuffled_c.nwe=0;
            for(size_t i=0; i<POLY_PART; i++)
            {
                auto nc = ids[i];

                // Make sure that the ids are still the same
                assert( sorted_c.contacts[nc].id == nc );

                auto age_part=sorted_c.contacts[nc].age;
                shuffled_c.ni[age_part]++;
                if(sorted_c.contacts[nc].weekend) shuffled_c.nwe++;

                shuffled_c.contacts[i]=sorted_c.contacts[nc];
            }

            return shuffled_c;
        }

        Eigen::MatrixXd to_symmetric_matrix( const contacts_t &c, 
                const data::age_data_t &age_data )
        {
            double ww[POLY_PART], mij[49], w_norm[7], cij[49], cij_pro;

            /*update of the weights*/
            for(size_t i=0;i<7;i++)
                w_norm[i]=0;

            for(size_t i=0;i<49;i++)
                mij[i]=0;

            for(size_t i=0; i<POLY_PART; i++)
            {
                int age_part=c.contacts[i].age;
                int AG_part=c.contacts[i].AG;
                if(!c.contacts[i].weekend)
                    ww[i]=(double)age_data.age_sizes[age_part]/c.ni[age_part]*5/(POLY_PART-c.nwe);
                else
                    ww[i]=(double)age_data.age_sizes[age_part]/c.ni[age_part]*2/c.nwe;

                w_norm[AG_part]+=ww[i];
                mij[7*AG_part]+=c.contacts[i].N1*ww[i];
                mij[7*AG_part+1]+=c.contacts[i].N2*ww[i];
                mij[7*AG_part+2]+=c.contacts[i].N3*ww[i];
                mij[7*AG_part+3]+=c.contacts[i].N4*ww[i];
                mij[7*AG_part+4]+=c.contacts[i].N5*ww[i];
                mij[7*AG_part+5]+=c.contacts[i].N6*ww[i];
                mij[7*AG_part+6]+=c.contacts[i].N7*ww[i];
            }

            /*Compute the contact matrix*/
            for(size_t i=0; i<49; i++)
            {
                if(w_norm[i/7]>0)
                    mij[i]/=w_norm[i/7];
                cij[i]=mij[i]/age_data.age_group_sizes[i%7];
            }

            Eigen::MatrixXd contact_regular( 7, 7 );
            for(size_t i=0; i<contact_regular.rows(); i++)
            {
                contact_regular(i,i)=cij[i*7+i];
                for(size_t j=0;j<i;j++)
                {
                    cij_pro=(cij[i*7+j]+cij[j*7+i])/2;
                    contact_regular(i,j)=cij_pro;
                    contact_regular(j,i)=cij_pro;
                }
            }

            return contact_regular;
        }
    };
};
