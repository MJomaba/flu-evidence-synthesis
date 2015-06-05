#ifndef FLU_DATA_HH
#define FLU_DATA_HH

#include<string>
#include<vector>
#include<array>

#include "state.hh"

namespace flu
{
    namespace data
    {
        std::vector<double> load_population( const std::array<size_t,7>
                &group_sizes );

        struct age_data_t
        {
            /// Population per year/age
            std::array<size_t,85> age_sizes;

            /// Population per age group
            std::array<size_t,7> age_group_sizes;
        };

        /// Load size of population per age and per age group
        age_data_t load_age_data( const std::string &path );
    };
};

#endif
