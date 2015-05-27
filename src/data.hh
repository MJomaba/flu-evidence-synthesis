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
    };
};

#endif
