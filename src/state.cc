#include "state.hh"

#include "io.hh"

namespace flu {

state_t load_state( const std::string &file_path )
{
    auto f_init = read_file( file_path );
    state_t state;
    return state;
}

};
