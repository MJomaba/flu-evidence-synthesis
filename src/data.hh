#ifndef FLU_DATA_HH
#define FLU_DATA_HH

#include<string>
#include<vector>

#include "state.hh"

namespace flu
{
    namespace data
    {
        std::vector<double> load_population( const std::string &filename );

        std::vector<double> load_contact_regular( const std::string& contacts_filepath, const state_t &state, int *age_sizes, int *AG_sizes  );
    };
};

#endif
