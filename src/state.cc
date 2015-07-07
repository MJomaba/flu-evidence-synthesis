#include<fstream>

#include "io.hh"
#include "state.hh"
#include "json.hh"

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
        //bbuild << "bootstrapped_contact_ids" << state.contact_ids;
        //
        bbuild << "likelihood" << state.likelihood;
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
        std::ifstream ifs( file_path );
        std::string content( (std::istreambuf_iterator<char>(ifs) ),
                (std::istreambuf_iterator<char>() ) );

        auto state = json::from_json_string<state_t>( content );

        return state;
    }

    void save_state_json( const state_t &state, 
            const std::string &file_path )
    {
        std::ofstream of( file_path );
        std::string json_str = json::to_json_string<state_t>( state ) ;
        of << json_str;
        of.close();
    }
};
