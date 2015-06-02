#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include <boost/numeric/ublas/io.hpp>

#include "model.hh"
#include "state.hh"
#include "data.hh"
#include "contacts.hh"
#include "vaccine.hh"
#include "json.hh"

using namespace flu;

int main(int argc, char *argv[])
{
    // Command line options
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;
   	po::options_description desc( "Analyses posterior distribution of the contact matrix\nUsage: contact_matrix --data-path [DIR]" );

    std::string data_path = "./";

    desc.add_options()
        ("help,h", "This message.")
        ("data-path,d", po::value<std::string>( &data_path ), "Path to samples/ directory)")
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

    auto age_data = data::load_age_data( data_path + "age_sizes.txt" );

    auto c = contacts::load_contacts( 
            data_path + "contacts_for_inference.txt" );


    auto contact_matrix = contacts::to_symmetric_matrix( 
            c, age_data );

    std::cout << "Original matrix" << std::endl;
    std::cout << contact_matrix << std::endl;

    for( auto & k : ks ) 
    {
        std::string kpadded = boost::lexical_cast<std::string>( k );
        while (kpadded.size() < 4)
            kpadded = "0" + kpadded;

        std::cout << kpadded << std::endl;

        auto state = load_state_json( data_path + "samples/z_hyper" 
                + kpadded
                + ".stm" );

        contact_matrix = contacts::to_symmetric_matrix( 
                contacts::shuffle_by_id( c, 
                    state.contact_ids ), age_data );

        std::cout << contact_matrix << std::endl;

/*        mongo::BSONEmitter bbuild;
        bbuild << "bla" << contact_matrix;
        auto bobj = bbuild.obj(); 
        auto json_str = bobj.jsonString( mongo::Strict, 1 );*/
 
        //std::string json_str = json::to_json_string( contact_matrix );
        //std::cout << json_str << std::endl;

    }

    return 0;
}
