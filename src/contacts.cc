#include "contacts.hh"

#include <stdio.h>
#include <cassert>

#include "io.hh"
#include "state.hh"

namespace flu
{
    namespace contacts {
        contacts_t load_contacts( const std::string &path )
        {
            auto contacts_PM = read_file( path );
            contacts_t c;

            for(size_t i=0; i<90; i++)
                c.ni[i]=0;

            c.nwe=0;

            /*Loading of the participants with their number of contacts from Polymod*/
            for(size_t i=0; i<POLY_PART; i++)
            {
                save_fscanf(contacts_PM,"%d %d %d %d %d %d %d %d %d", &c.contacts[i].age, &c.contacts[i].we, &c.contacts[i].N1, &c.contacts[i].N2, &c.contacts[i].N3, &c.contacts[i].N4, &c.contacts[i].N5, &c.contacts[i].N6, &c.contacts[i].N7);
                auto age_part=c.contacts[i].age;
                c.ni[age_part]++;
                if(c.contacts[i].we>0) c.nwe++;
                c.contacts[i].AG=0;
                if(age_part>0) c.contacts[i].AG++;
                if(age_part>4) c.contacts[i].AG++;
                if(age_part>14) c.contacts[i].AG++;
                if(age_part>24) c.contacts[i].AG++;
                if(age_part>44) c.contacts[i].AG++;
                if(age_part>64) c.contacts[i].AG++;
                c.contacts[i].id = i;
            }

            return c;
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
                if(sorted_c.contacts[nc].we>0) shuffled_c.nwe++;

                shuffled_c.contacts[i]=sorted_c.contacts[nc];
            }

            return shuffled_c;
        }

        std::vector<double> to_symmetric_matrix( const contacts_t &c, int *age_sizes, int *AG_sizes )
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
                if(c.contacts[i].we==0)
                    ww[i]=(double)age_sizes[age_part]/c.ni[age_part]*5/(POLY_PART-c.nwe);
                else
                    ww[i]=(double)age_sizes[age_part]/c.ni[age_part]*2/c.nwe;

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

            return contact_regular;
        }
    };
};
