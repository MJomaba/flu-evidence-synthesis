#ifndef FLU_DATA_HH
#define FLU_DATA_HH

#include<string>
#include<vector>
#include<array>

#include "state.h"

namespace flu
{
    namespace data
    {
        /** 
         * \brief Separate the groups sizes into risk groups
         *
         * Groups are: low risk, high risk and pregnant women
         */
        std::vector<double>  separate_into_risk_groups( 
                const std::vector<size_t> &group_sizes );

        /**
         * \brief Group population size according to age groups
         */
        std::vector<size_t> group_age_data( 
                const std::vector<size_t> &age_sizes, 
                const std::vector<size_t> &limits );

        struct age_data_t
        {
            /// Population per year/age
            std::vector<size_t> age_sizes;

            /// Population per age group
            std::vector<size_t> age_group_sizes;
        };
    }
}

#endif
