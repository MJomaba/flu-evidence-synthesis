#include "state.hh"

#include "io.hh"

namespace flu {

state_t load_state( const std::string &file_path, 
        const size_t number_age_groups, const size_t dim_poly_part )
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

    save_fgets(sbuffer, 100, f_init);
    for(size_t i=0;i<52;i++)
    {
        save_fgets(sbuffer, 150, f_init);
        sscanf(sbuffer,"%lf %lf %lf %lf %lf",&state.positivity_ij[i*5],&state.positivity_ij[i*5+1],&state.positivity_ij[i*5+2],&state.positivity_ij[i*5+3],&state.positivity_ij[i*5+4]);
    }

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
    sscanf(sbuffer,"%lf %lf %lf %lf %lf", &state.parameters.epsilon[0],&state.parameters.epsilon[1],&state.parameters.epsilon[2],&state.parameters.epsilon[3],&state.parameters.epsilon[4]);

    save_fgets(sbuffer, 100, f_init);
    save_fgets(sbuffer, 100, f_init);
	sscanf(sbuffer,"%lf", &state.parameters.psi);

    save_fgets(sbuffer, 100, f_init);

    for(size_t i=0; i<dim_poly_part; i++)
    {
        int nc;
        save_fgets(sbuffer, 10, f_init);
        sscanf(sbuffer,"%d",&nc);
        state.number_contacts.push_back( nc );
    }

    return state;
}

};
