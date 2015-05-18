#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

#include "model.hh"
#include "state.hh"
#include "data.hh"

using namespace flu;

int main(int argc, char *argv[])
{
    // Command line options
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;
   	po::options_description desc( "Usage: inference --data-path [DIR]" );

    std::string data_path = "./";

    desc.add_options()
        ("help,h", "This message.")
        ("data-path,d", po::value<std::string>( &data_path ), "Path to the posterior.txt file (and samples/ directory)")
        ;

    po::variables_map vm;
    po::store( 
            po::command_line_parser( argc, argv ).options(desc).run(),
            vm );

    try {
        po::notify( vm );
    } catch (po::required_option e) {
        std::cout << e.what() << std::endl << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    } 

    data_path = (fs::canonical(
                fs::complete( data_path ) )).native() + "/";

    fs::directory_iterator it(data_path + "samples/");
    fs::directory_iterator endit;

    std::vector<size_t> ks;

    static const boost::regex kRegex("z_hyper(\\d+)\\.stm");
    while(it != endit)
    {
        if(fs::is_regular_file(*it) 
          && it->path().extension() == ".stm" ) 
        {
            boost::match_results<std::string::const_iterator> results;
            if (boost::regex_match( 
                        it->path().filename().native(), results, kRegex))
                ks.push_back( boost::lexical_cast<size_t>( results[1] ) );
        }
        ++it;
    }

    sort( ks.begin(), ks.end() );



    FILE *Scen1FS, *Scen2FS;
    Scen1FS=write_file(data_path + "scenarii/Scenario_vaccination_final_size.txt");

    /*opens the 2nd scenarioFile*/
    Scen2FS=write_file(data_path + "scenarii/Scenario_no_vaccination_final_size.txt");

    /*load the size of the age groups for the model that year*/
    auto pop_vec = data::load_population( data_path 
            + "age_groups_model.txt" );

    //double current_contact_regular[NAG2];
    
    //int n_scenarii;

    for( auto & k : ks ) 
    {
        auto state = load_state( data_path + "samples/z_hyper" 
                + boost::lexical_cast<std::string>( k )
                + ".stm", NAG, POLY_PART );

        /*translate into an initial infected population*/
        double init_inf[NAG];
        for(size_t i=0;i<NAG;i++)
            init_inf[i]=pow(10,state.parameters.init_pop);

        //save_scenarii(Scen1FS, Scen2FS, pop_vec, init_inf, state, current_contact_regular, n_scenarii, tab_cal, tab_VE, data_path, &First_write);
    }
    fclose(Scen1FS);
    fclose(Scen2FS);

    return 0;
}
