#ifndef FLU_CONTACTS_HH
#define FLU_CONTACTS_HH

#include<array>
#include<vector>
#include<string>

#include "state.h"
#include "data.h"

#include<Eigen/Core>

namespace flu {

    namespace contacts {

        /// Create a contact, holding its age, and contacts (we, N1, N2, etc)..
        struct contact_t
        {
            size_t id;
            int age, AG;

            std::vector<int> N;

            /// Is this data from the weekend
            bool weekend;
        };

        /// Struct to keep all contacts and some metadata
        struct contacts_t
        {
            std::vector<contact_t> contacts;

            /// Total number of participants of each age group
            int ni[90];

            /// Total amount of contact data for weekend
            size_t nwe;
        };

        /// Convert polymod table (matrix) to contacts_t struct
        contacts_t table_to_contacts(
                const Eigen::MatrixXi &conMatrix,
                const std::vector<size_t> &limits );

        /// Bootstrap the given contacts, by shuffling back the given no of contacts from the original data
        contacts_t bootstrap_contacts( contacts_t&& bootstrap,
                const contacts_t &original,
                size_t no ); 

         /**
         * \brief Shuffle given contacts according to id. Assumes the given
         * contacts are already sorted
         */
        contacts_t shuffle_by_id( const contacts_t &sorted_c, const std::vector<size_t> &ids );

        Eigen::MatrixXd to_symmetric_matrix( 
                const contacts_t &contacts, 
                const data::age_data_t &age_data );
    }

}
#endif

