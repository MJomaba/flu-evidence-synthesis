#include "state.hh"

#include "io.hh"

namespace flu {

State load_state( const std::string &file_path )
{
    auto f_init = read_file( file_path );
    State state;
    return state;
}

};
