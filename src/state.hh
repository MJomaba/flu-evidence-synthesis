#ifndef STATE_HH 
#define STATE_HH

#include<string>

/// Keeps the current state of the model/mcmc
struct State 
{
};

/// Load an (initial) state from a file
State load_state( const std::string &file_path );

#endif
