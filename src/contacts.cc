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
                auto alea1=(size_t) R::runif(0,bootstrap.contacts.size());
                auto alea2=(size_t) R::runif(0,bootstrap.contacts.size());

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
            for(size_t i=0; i<sorted_c.contacts.size(); i++)
            {
                auto nc = ids[i];

                // Make sure that the ids are still the same
                assert( sorted_c.contacts[nc].id == nc );

                auto age_part=sorted_c.contacts[nc].age;
                shuffled_c.ni[age_part]++;
                if(sorted_c.contacts[nc].weekend) shuffled_c.nwe++;

                shuffled_c.contacts.push_back(sorted_c.contacts[nc]);
            }

            return shuffled_c;
        }

        Eigen::MatrixXd to_symmetric_matrix( const contacts_t &c, 
                const data::age_data_t &age_data )
        {
            Eigen::VectorXd ww = Eigen::VectorXd( c.contacts.size() );
            auto nag = age_data.age_group_sizes.size();
            Eigen::MatrixXd mij = Eigen::MatrixXd::Zero(nag, nag);
            Eigen::VectorXd w_norm = Eigen::VectorXd::Zero(nag);
            Eigen::MatrixXd cij(nag, nag);

            for(size_t i=0; i<c.contacts.size(); i++)
            {
                int age_part=c.contacts[i].age;
                int AG_part=c.contacts[i].AG;
                if(!c.contacts[i].weekend)
                    ww[i]=(double)age_data.age_sizes[age_part]/c.ni[age_part]*5/(c.contacts.size()-c.nwe);
                else
                    ww[i]=(double)age_data.age_sizes[age_part]/c.ni[age_part]*2/c.nwe;
                w_norm[AG_part]+=ww[i];
                for (size_t j=0; j<c.contacts[i].N.size(); ++j)
                    mij(AG_part,j)+=c.contacts[i].N[j]*ww[i];
            }

            /*Compute the contact matrix*/
            for(int i=0; i<mij.cols(); i++)
            {
                for(int j=0; j<mij.rows(); j++)
                {
                    if(w_norm[j]>0)
                        mij(j,i)/=w_norm[j];
                    cij(j,i)=mij(j,i)/age_data.age_group_sizes[i];
                }

            }

            Eigen::MatrixXd contact_regular( nag, nag );
            for(int i=0; i<contact_regular.rows(); i++)
            {
                contact_regular(i,i)=cij(i,i);
                for(int j=0;j<i;j++)
                {
                    auto cij_pro=(cij(i,j)+cij(j,i))/2;
                    contact_regular(i,j)=cij_pro;
                    contact_regular(j,i)=cij_pro;
                }
            }

            return contact_regular;
        }
    }
}
