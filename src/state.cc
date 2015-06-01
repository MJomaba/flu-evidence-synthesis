#include<fstream>

#include "client/redef_macros.h"
#include "db/json.h"

#include "io.hh"
#include "state.hh"

gsl_rng * r;

namespace flu {

    mongo::BSONEmitter &operator<<(
            mongo::BSONEmitter &bbuild, const parameter_set &pars )
    {
        bbuild << "epsilon" << pars.epsilon;
        bbuild << "psi" << pars.psi;
        bbuild << "transmissibility" << pars.transmissibility;
        bbuild << "susceptibility" << pars.susceptibility;
        bbuild << "init_pop" << pars.init_pop;
        return bbuild;
    }

    void operator>>( const mongo::BSONElement &el,
            parameter_set &pars )
    {
        el["epsilon"] >> pars.epsilon;
        el["psi"] >> pars.psi;
        el["transmissibility"] >> pars.transmissibility;
        el["susceptibility"] >> pars.susceptibility;
        el["init_pop"] >> pars.init_pop;
    }

    mongo::BSONEmitter &operator<<(
            mongo::BSONEmitter &bbuild, const state_t &state )
    {
        bbuild << "parameters" << state.parameters;
        bbuild << "time_infectious" << state.time_infectious;
        bbuild << "time_latent" << state.time_latent;
        bbuild << "positivity" << state.positivity_ij;
        std::vector<int> conv;
        conv.resize( state.contact_ids.size() );

        // For readability of the json output, convert to int first
        std::transform( state.contact_ids.begin(),
                state.contact_ids.end(), conv.begin(),
                [](size_t i) {return (int)i;} );
        bbuild << "bootstrapped_contact_ids" << conv;
        return bbuild;
    }

    void operator>>( const mongo::BSONElement &el,
            state_t &state )
    {
        el["parameters"] >> state.parameters;
        el["time_infectious"] >> state.time_infectious;
        el["time_latent"] >> state.time_latent;
        el["positivity"] >> state.positivity_ij;
        el["bootstrapped_contact_ids"] >> state.contact_ids;
    }


    state_t load_state_json( const std::string &file_path ) 
    {
        state_t state;

        std::ifstream ifs( file_path );
        std::string content( (std::istreambuf_iterator<char>(ifs) ),
                (std::istreambuf_iterator<char>() ) );

        mongo::BSONObj json_state = mongo::fromjson( content );

        json_state >> state;
        return state;
    }

    void save_state_json( const state_t &state, 
            const std::string &file_path )
    {
        mongo::BSONEmitter bbuild;
        bbuild << state;

        std::ofstream of( file_path );
        of << bbuild.obj().jsonString( mongo::Strict, 1 );
        of.close();
    }

};
