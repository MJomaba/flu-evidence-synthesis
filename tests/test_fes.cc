#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <iostream>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "catch.hpp"

#include "model.hh"
#include "trace.hh"

namespace io = boost::iostreams;
namespace fs = boost::filesystem;

bool equal_file_contents( const std::string &path, 
        const std::string &other_path )
{
    io::mapped_file_source f1( path );
    io::mapped_file_source f2( other_path );

    if ( f1.size() > 0 
            && f1.size() == f2.size() 
            && std::equal( f1.data(), 
                f1.data() + f1.size(), f2.data() ) 
       ) 
        return true;
    else
        return false;
}

TEST_CASE( "Run a short test run", "[full]" ) 
{
    // Create array with all names to be watched
    std::vector<std::string> files;
    files.push_back( "posterior.txt" );

    // Add sample files to check
    for( size_t i = 1010; i < 2001; i+= 10 )
    {
        files.push_back( "samples/z_hyper" 
                + boost::lexical_cast<std::string>( i ) + ".stm" );
    }

    // Setup the scenarii for watching
    std::string pre = "scenarii/Scenario_";
    std::string post = "_final_size.txt";
    files.push_back( pre + "vaccination" + post );
    files.push_back( pre + "no_vaccination" + post );
    for ( size_t i = 0; i<35; ++i )
        files.push_back( pre + boost::lexical_cast<std::string>( i )
                + post );

    // Start with clean slate
    for( auto &fname : files )
    {
        fs::remove( "../data/" + fname );
    }


    REQUIRE( system( "bin/flu-evidence-synthesis -d ../data/ --burn-in 1000 --chain-length 1000 --thinning 1" ) == 0 );

    for( auto &fname : files )
    {
        REQUIRE( fs::exists("../data/" + fname) );

        auto full_path = "../data/" + fname;
        auto test_path = "./tests/test_data/" + fname;

        REQUIRE( equal_file_contents( 
                    full_path, test_path ) );
    }
}

TEST_CASE( "Run inference", "[hide]" )
{
    // Create array with all names to be watched
    std::vector<std::string> files;

    // Setup the scenarii for watching
    std::string pre = "scenarii/Scenario_";
    std::string post = "_final_size.txt";
    files.push_back( pre + "vaccination" + post );
    files.push_back( pre + "no_vaccination" + post );
    for ( size_t i = 0; i<35; ++i )
        files.push_back( pre + boost::lexical_cast<std::string>( i )
                + post );

    // Start with clean slate
    for( auto &fname : files )
    {
        fs::remove( "../data/" + fname );
    }


    REQUIRE( system( "bin/inference -d ../data/" ) == 0 );

    for( auto &fname : files )
    {
        REQUIRE( fs::exists("../data/" + fname) );

        auto full_path = "../data/" + fname;
        auto test_path = "./tests/test_data/" + fname;

        REQUIRE( equal_file_contents( 
                    full_path, test_path ) );
    }
}

TEST_CASE( "Posterior file should be readable as a trace" "[trace]" )
{
    std::string fName = "./tests/test_data/posterior.txt";

    std::ifstream infile;
    infile.open( fName );

    std::vector<std::vector<double> > samples 
        = flu::trace::read_trace_per_sample( infile, 1 );

    REQUIRE( samples.size() == 1000 );
    REQUIRE( samples[0].size() == 12 );

    infile.close();
}
