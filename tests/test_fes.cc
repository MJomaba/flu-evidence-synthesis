#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <iostream>
#include <boost/iostreams/device/mapped_file.hpp>

#include "catch.hpp"

#include "model.hh"

namespace io = boost::iostreams;

TEST_CASE( "Run a short test run", "[full]" ) {
    REQUIRE( system( "bin/flu-evidence-synthesis -d ../data/ --burn-in 1000 --chain-length 1000 --thinning 1" ) == 0 );
    io::mapped_file_source f1( "../data/posterior.txt" );
    io::mapped_file_source f2( "./tests/short_posterior_sample.txt" );

    REQUIRE( f1.size() > 0 );
    REQUIRE( f1.size() == f2.size() );
    REQUIRE( std::equal( f1.data(), f1.data() + f1.size(), f2.data() ) );
}
