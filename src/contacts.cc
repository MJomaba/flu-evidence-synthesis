#include "contacts.hh"

#include <stdio.h>

#include "io.hh"
#include "state.hh"

namespace flu
{
    namespace contacts {
        std::vector<double> load_contact_regular( const std::string& contacts_filepath, const state_t &state, int *age_sizes, int *AG_sizes, const contacts_t &c )
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
