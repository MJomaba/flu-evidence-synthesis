#include "state.hh"

#include "io.hh"

namespace flu {

state_t load_state( const std::string &file_path, 
        const size_t number_age_groups )
{
    char sbuffer[300];
    auto f_init = read_file( file_path );

    state_t state;

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf", &state.time_latent, &state.time_infectious);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf", &state.parameters.init_pop);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf", &state.parameters.transmissibility);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf %lf %lf", &state.parameters.susceptibility[0],&state.parameters.susceptibility[1],&state.parameters.susceptibility[2],&state.parameters.susceptibility[3],&state.parameters.susceptibility[4],&state.parameters.susceptibility[5],&state.parameters.susceptibility[6]);

    return state;
}

};
