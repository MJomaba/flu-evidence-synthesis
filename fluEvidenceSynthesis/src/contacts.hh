#ifndef FLU_CONTACTS_HH
#define FLU_CONTACTS_HH

#include<array>
#include<vector>
#include<string>

#include <boost/numeric/ublas/matrix.hpp>

#include "state.hh"
#include "data.hh"

namespace flu {

    namespace contacts {
        namespace bu = boost::numeric::ublas;

        /// Create a contact, holding its age, and contacts (we, N1, N2, etc)..
        struct contact_t
        {
            size_t id;
            int age, N1, N2, N3, N4, N5, N6, N7, AG;

            /// Is this data from the weekend
            bool weekend;
        };

        /// Struct to keep all contacts and some metadata
        struct contacts_t
        {
            std::array<contact_t,POLY_PART> contacts;

            /// Total number of participants of each age group
            int ni[90];

            /// Total amount of contact data for weekend
            size_t nwe;
        };

        /// Bootstrap the given contacts, by shuffling back the given no of contacts from the original data
        contacts_t bootstrap_contacts( contacts_t&& bootstrap,
                const contacts_t &original,
                size_t no ); 

         /**
         * \brief Shuffle given contacts according to id. Assumes the given
         * contacts are already sorted
         */
        contacts_t shuffle_by_id( const contacts_t &sorted_c, const std::vector<size_t> &ids );

        bu::matrix<double> to_symmetric_matrix( 
                const contacts_t &contacts, 
                const data::age_data_t &age_data );
    };

};
#endif

