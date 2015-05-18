#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

#include "model.hh"
#include "state.hh"

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
    FILE *f_pop_model;
    double pop_vec[21];
    save_fscanf(f_pop_model,"%lf %lf %lf %lf %lf %lf %lf",&pop_vec[0],&pop_vec[1],&pop_vec[2],&pop_vec[3],&pop_vec[4],&pop_vec[5],&pop_vec[6]);
    /*high risk*/
    pop_vec[7]=pop_vec[0]*0.021; /*2.1% in the <1 */
    pop_vec[8]=pop_vec[1]*0.055; /*5.5% in the 1-4 */
    pop_vec[9]=pop_vec[2]*0.098; /*9.8% in the 5-14 */
    pop_vec[10]=pop_vec[3]*0.087; /*8.7% in the 15-24 */
    pop_vec[11]=pop_vec[4]*0.092; /*9.2% in the 25-44 */
    pop_vec[12]=pop_vec[5]*0.183; /*18.3% in the 45-64 */
    pop_vec[13]=pop_vec[6]*0.45; /*45% in the 65+ */

    /*pregnant women (pregnant women not considered in this model*/
    pop_vec[14]=0;
    pop_vec[15]=0;
    pop_vec[16]=0;
    pop_vec[17]=0;
    pop_vec[18]=0;
    pop_vec[19]=0;
    pop_vec[20]=0;

    /*low risk*/
    pop_vec[0]-=pop_vec[7];
    pop_vec[1]-=pop_vec[8];
    pop_vec[2]-=pop_vec[9];
    pop_vec[3]-=pop_vec[10]+pop_vec[17];
    pop_vec[4]-=pop_vec[11]+pop_vec[18];
    pop_vec[5]-=pop_vec[12];
    pop_vec[6]-=pop_vec[13];
    fclose( f_pop_model );

    double curr_init_inf[NAG];
    //double current_contact_regular[NAG2];
    
    //int n_scenarii;

    for( auto & k : ks ) 
    {
        auto state = load_state( data_path + "samples/z_hyper" 
                + boost::lexical_cast<std::string>( k )
                + ".stm", NAG, POLY_PART );

        /*translate into an initial infected population*/
        for(size_t i=0;i<NAG;i++)
            curr_init_inf[i]=pow(10,state.parameters.init_pop);

        //save_scenarii(Scen1FS, Scen2FS, pop_vec, curr_init_inf, state, current_contact_regular, n_scenarii, tab_cal, tab_VE, data_path, &First_write);
        std::cout << k << std::endl;
    }
    fclose(Scen1FS);
    fclose(Scen2FS);

    return 0;
}
